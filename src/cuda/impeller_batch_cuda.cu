// SPDX-License-Identifier: GPL-3.0-or-later
#include <turbocompgpu/gpu_impeller_batch_cuda.hpp>

#include <friction/colebrook_cuda.cuh>
#include <cuda_runtime.h>

#include <cmath>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>

namespace turbocompgpu::cuda {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

struct Residual {
    double f0, f1, norm;
    double p, rho, vm, whirl, work, torque, er, sr;
    double friction_factor, friction_entropy;
    thermomodel::detail::PRMixtureCudaState thermo;
    bool valid;
};

__device__ bool valid_input(const ImpellerBatchInput& x,
                            const ImpellerBatchOptions& options)
{
    return isfinite(x.inlet_temperature_K) && x.inlet_temperature_K > 0.0 &&
        isfinite(x.inlet_pressure_Pa) && x.inlet_pressure_Pa > 0.0 &&
        isfinite(x.mass_flow_kg_s) && x.mass_flow_kg_s > 0.0 &&
        isfinite(x.speed_rpm) && x.speed_rpm > 0.0 &&
        isfinite(x.slip_factor) && x.slip_factor > 0.0 && x.slip_factor <= 1.0 &&
        isfinite(x.inlet_prewhirl_m_s) &&
        isfinite(x.base_entropy_generation_J_kgK) && x.base_entropy_generation_J_kgK >= 0.0 &&
        isfinite(x.inlet_hub_diameter_m) && x.inlet_hub_diameter_m >= 0.0 &&
        isfinite(x.inlet_shroud_diameter_m) && x.inlet_shroud_diameter_m > x.inlet_hub_diameter_m &&
        isfinite(x.inlet_blockage) && x.inlet_blockage >= 0.0 && x.inlet_blockage < 1.0 &&
        isfinite(x.exit_diameter_m) && x.exit_diameter_m > 0.0 &&
        isfinite(x.exit_width_m) && x.exit_width_m > 0.0 &&
        isfinite(x.exit_blockage) && x.exit_blockage >= 0.0 && x.exit_blockage < 1.0 &&
        isfinite(x.inlet_impeller_diameter_m) && x.inlet_impeller_diameter_m > 0.0 &&
        x.inlet_impeller_diameter_m < x.exit_diameter_m &&
        isfinite(x.exit_blade_angle_deg) && x.exit_blade_angle_deg > 0.0 && x.exit_blade_angle_deg <= 90.0 &&
        isfinite(x.specific_work_scale_J_kg) && x.specific_work_scale_J_kg > 0.0 &&
        isfinite(x.specific_entropy_scale_J_kgK) && x.specific_entropy_scale_J_kgK > 0.0 &&
        isfinite(x.pressure_scale_Pa) && x.pressure_scale_Pa > 0.0 &&
        isfinite(x.passage_length_m) && x.passage_length_m >= 0.0 &&
        isfinite(x.hydraulic_diameter_m) && x.hydraulic_diameter_m > 0.0 &&
        (options.use_iso20765_viscosity ||
         (isfinite(x.dynamic_viscosity_Pa_s) && x.dynamic_viscosity_Pa_s > 0.0)) &&
        isfinite(x.relative_roughness) && x.relative_roughness >= 0.0 && x.relative_roughness < 3.7 &&
        isfinite(x.initial_temperature_K) && x.initial_temperature_K > 0.0 &&
        isfinite(x.initial_log_pressure);
}

__device__ Residual evaluate_scaled(const ImpellerBatchInput& x,
    const ImpellerBatchOptions& o,
    const thermomodel::detail::PRMixtureCudaParameters& mix,
    const transport::Iso20765Parameters& viscosity_model,
    const thermomodel::detail::PRMixtureCudaState& inlet,
    double y0, double y1)
{
    Residual r{};
    const double T = y0 * o.temperature_scale_K;
    const double p = x.pressure_scale_Pa * exp(y1);
    if (!(isfinite(T) && T > 0.0 && isfinite(p) && p > 0.0)) return r;
    r.thermo = thermomodel::detail::pr_mixture_evaluate_Tp_cuda(mix, T, p);
    if (r.thermo.status != thermomodel::detail::PRMixtureCudaStatus::Success ||
        !(isfinite(r.thermo.rho) && r.thermo.rho > 0.0)) return r;

    const double inlet_area = (pi/4.0) *
        (x.inlet_shroud_diameter_m*x.inlet_shroud_diameter_m -
         x.inlet_hub_diameter_m*x.inlet_hub_diameter_m) * (1.0-x.inlet_blockage);
    const double exit_area = pi*x.exit_diameter_m*x.exit_width_m*(1.0-x.exit_blockage);
    const double vm1 = x.mass_flow_kg_s/(inlet.rho*inlet_area);
    r.vm = x.mass_flow_kg_s/(r.thermo.rho*exit_area);
    const double omega = 2.0*pi*x.speed_rpm/60.0;
    const double U1 = omega*x.inlet_impeller_diameter_m/2.0;
    const double U2 = omega*x.exit_diameter_m/2.0;
    const double beta = x.exit_blade_angle_deg*pi/180.0;
    r.whirl = x.slip_factor*(U2-r.vm*cos(beta)/sin(beta));
    r.work = U2*r.whirl-U1*x.inlet_prewhirl_m_s;
    const double h01 = inlet.h+0.5*(vm1*vm1+x.inlet_prewhirl_m_s*x.inlet_prewhirl_m_s);
    const double h02 = r.thermo.h+0.5*(r.vm*r.vm+r.whirl*r.whirl);
    r.er = (h02-h01)-r.work;
    r.torque = x.mass_flow_kg_s*(x.exit_diameter_m*r.whirl/2.0 -
        x.inlet_impeller_diameter_m*x.inlet_prewhirl_m_s/2.0);

    const double relative_speed = hypot(vm1, x.inlet_prewhirl_m_s);
    const double mu_inlet = o.use_iso20765_viscosity
        ? transport::iso20765_5_dynamic_viscosity_Pa_s(
              viscosity_model, inlet.T, inlet.p, inlet.rho)
        : x.dynamic_viscosity_Pa_s;
    if (!(isfinite(mu_inlet) && mu_inlet > 0.0)) return r;
    const double Re = inlet.rho*relative_speed*x.hydraulic_diameter_m/mu_inlet;
    const auto friction = friction::cuda::colebrook_turbulent(Re, x.relative_roughness);
    if (friction.status != friction::cuda::ColebrookStatus::Success) return r;
    r.friction_factor = friction.darcy_friction_factor;
    const double friction_dh = friction.darcy_friction_factor *
        (x.passage_length_m/x.hydraulic_diameter_m) * 0.5*relative_speed*relative_speed;
    r.friction_entropy = friction_dh/x.inlet_temperature_K;
    const double entropy_generation = x.base_entropy_generation_J_kgK+r.friction_entropy;
    r.sr = (r.thermo.s-inlet.s)-entropy_generation;
    r.f0 = r.er/x.specific_work_scale_J_kg;
    r.f1 = r.sr/x.specific_entropy_scale_J_kgK;
    r.norm = hypot(r.f0,r.f1);
    r.p = p;
    r.rho = r.thermo.rho;
    r.valid = isfinite(r.norm) && isfinite(r.torque) && isfinite(r.work);
    return r;
}

struct DiffuserRhs {
    double dp, dT, dv, ds, dL;
    bool valid;
};

struct DiffuserResult {
    double T, p, rho, vm, whirl;
    double energy_residual, entropy_change, entropy_integral;
    double entropy_residual, angular_residual, wall_torque_integral;
    double mass_residual;
    bool valid;
};

__device__ DiffuserRhs diffuser_rhs(
    const thermomodel::detail::PRMixtureCudaParameters& mix,
    const transport::Iso20765Parameters& viscosity_model,
    bool use_iso20765_viscosity,
    double radius, double width, double mass,
    double fixed_dynamic_viscosity, double relative_roughness,
    double p, double T, double tangential_velocity)
{
    DiffuserRhs out{};
    const auto state=thermomodel::detail::pr_mixture_evaluate_Tp_cuda(mix,T,p);
    if (state.status!=thermomodel::detail::PRMixtureCudaStatus::Success ||
        !(isfinite(state.rho) && state.rho>0.0 && isfinite(state.h))) return out;

    const double area=2.0*pi*radius*width;
    const double u=mass/(state.rho*area);
    if (!(isfinite(u) && u>0.0)) return out;
    const double speed=hypot(u,tangential_velocity);
    // Evaluate local transport viscosity, Re, and Colebrook at each RK stage.
    const double hydraulic_diameter=2.0*width;
    const double dynamic_viscosity = use_iso20765_viscosity
        ? transport::iso20765_5_dynamic_viscosity_Pa_s(
              viscosity_model,T,p,state.rho)
        : fixed_dynamic_viscosity;
    if (!(isfinite(dynamic_viscosity) && dynamic_viscosity > 0.0)) return out;
    const double reynolds=state.rho*speed*hydraulic_diameter/dynamic_viscosity;
    const auto friction_result=friction::cuda::colebrook_turbulent(
        reynolds,relative_roughness);
    if (friction_result.status!=friction::cuda::ColebrookStatus::Success) return out;
    const double friction_factor=friction_result.darcy_friction_factor;
    const double friction_accel=friction_factor*speed/(4.0*width);

    // Port of the CPU analytical diffuser RHS. ThermoGPU's PR-mixture state
    // carries analytic rho and entropy derivatives; Gibbs gives h_p and h_T.
    const double rho_p=state.drho_dp_T;
    const double rho_T=state.drho_dT_p;
    const double h_p=1.0/state.rho+T*state.ds_dp_T;
    const double h_T=T*state.ds_dT_p;
    const double dv=-tangential_velocity/radius-
        friction_accel*tangential_velocity/u;
    const double area_gradient=1.0/radius; // constant-width passage
    const double B1=u*u*area_gradient+
        tangential_velocity*tangential_velocity/radius-friction_accel*u;
    const double B2=u*u*area_gradient-tangential_velocity*dv;
    const double A11=(1.0-u*u*rho_p)/state.rho;
    const double A12=-u*u*rho_T/state.rho;
    const double A21=h_p-u*u*rho_p/state.rho;
    const double A22=h_T-u*u*rho_T/state.rho;
    const double determinant=A11*A22-A12*A21;
    const double scale=fabs(A11*A22)+fabs(A12*A21);
    if (!(isfinite(determinant) && isfinite(scale) && scale>0.0 &&
          fabs(determinant)>1.0e-10*scale)) return out;

    out.dp=(B1*A22-A12*B2)/determinant;
    out.dT=(A11*B2-B1*A21)/determinant;
    out.dv=dv;
    out.ds=friction_factor*speed*speed*speed/(4.0*width*u*T);
    out.dL=-friction_accel*radius*tangential_velocity/u;
    out.valid=isfinite(out.dp)&&isfinite(out.dT)&&isfinite(out.dv)&&
        isfinite(out.ds)&&isfinite(out.dL);
    return out;
}

__device__ DiffuserResult integrate_diffuser(
    const ImpellerBatchInput& x,
    const thermomodel::detail::PRMixtureCudaParameters& mix,
    const transport::Iso20765Parameters& viscosity_model,
    bool use_iso20765_viscosity,
    const Residual& impeller)
{
    DiffuserResult out{};
    const double r2=x.exit_diameter_m/2.0;
    const double r3=x.diffuser_outlet_radius_m;
    const double width=x.exit_width_m;
    const int steps=x.diffuser_radial_steps;
    const double mass=x.mass_flow_kg_s;
    if (!(r3>r2 && width>0.0 && steps>0 && x.exit_blockage==0.0)) return out;

    const double inlet_area=2.0*pi*r2*width;
    const double expected_mass=impeller.thermo.rho*inlet_area*impeller.vm;
    if (fabs(expected_mass-mass)>1.0e-8*fmax(1.0,mass)) return out;
    const double h02=impeller.thermo.h+0.5*(impeller.vm*impeller.vm+
        impeller.whirl*impeller.whirl);

    double p=impeller.p, T=impeller.thermo.T, v=impeller.whirl;
    double entropy_integral=0.0, angular_integral=0.0;
    const double dr=(r3-r2)/static_cast<double>(steps);
    for (int i=0;i<steps;++i) {
        const double r=r2+static_cast<double>(i)*dr;
        const DiffuserRhs k1=diffuser_rhs(mix,viscosity_model,use_iso20765_viscosity,
            r,width,mass,x.dynamic_viscosity_Pa_s,x.relative_roughness,p,T,v);
        if (!k1.valid) return out;
        const DiffuserRhs k2=diffuser_rhs(mix,viscosity_model,use_iso20765_viscosity,
            r+0.5*dr,width,mass,x.dynamic_viscosity_Pa_s,x.relative_roughness,
            p+0.5*dr*k1.dp,T+0.5*dr*k1.dT,v+0.5*dr*k1.dv);
        if (!k2.valid) return out;
        const DiffuserRhs k3=diffuser_rhs(mix,viscosity_model,use_iso20765_viscosity,
            r+0.5*dr,width,mass,x.dynamic_viscosity_Pa_s,x.relative_roughness,
            p+0.5*dr*k2.dp,T+0.5*dr*k2.dT,v+0.5*dr*k2.dv);
        if (!k3.valid) return out;
        const DiffuserRhs k4=diffuser_rhs(mix,viscosity_model,use_iso20765_viscosity,
            r+dr,width,mass,x.dynamic_viscosity_Pa_s,x.relative_roughness,
            p+dr*k3.dp,T+dr*k3.dT,v+dr*k3.dv);
        if (!k4.valid) return out;
        p+=dr/6.0*(k1.dp+2.0*k2.dp+2.0*k3.dp+k4.dp);
        T+=dr/6.0*(k1.dT+2.0*k2.dT+2.0*k3.dT+k4.dT);
        v+=dr/6.0*(k1.dv+2.0*k2.dv+2.0*k3.dv+k4.dv);
        entropy_integral+=dr/6.0*(k1.ds+2.0*k2.ds+2.0*k3.ds+k4.ds);
        angular_integral+=dr/6.0*(k1.dL+2.0*k2.dL+2.0*k3.dL+k4.dL);
    }

    const auto state3=thermomodel::detail::pr_mixture_evaluate_Tp_cuda(mix,T,p);
    if (state3.status!=thermomodel::detail::PRMixtureCudaStatus::Success ||
        !(state3.rho>0.0 && isfinite(state3.rho))) return out;
    const double area3=2.0*pi*r3*width;
    const double vm3=mass/(state3.rho*area3);
    const double h03=state3.h+0.5*(vm3*vm3+v*v);
    const double entropy_change=state3.s-impeller.thermo.s;
    const double angular_change=r3*v-r2*impeller.whirl;
    out={T,p,state3.rho,vm3,v,h03-h02,entropy_change,entropy_integral,
         entropy_change-entropy_integral,angular_change-angular_integral,
         angular_integral,state3.rho*area3*vm3-mass,true};
    return out;
}

__device__ ImpellerBatchOutput failure(ImpellerBatchStatus status, int iteration)
{
    ImpellerBatchOutput out{};
    out.status = status;
    out.iterations = iteration;
    out.residual_norm = INFINITY;
    return out;
}

__global__ void solve_kernel(const ImpellerBatchInput* inputs,
    ImpellerBatchOutput* outputs, std::size_t n,
    thermomodel::detail::PRMixtureCudaParameters mixture,
    transport::Iso20765Parameters viscosity_model,
    ImpellerBatchOptions options)
{
    const std::size_t i = static_cast<std::size_t>(blockIdx.x)*blockDim.x+threadIdx.x;
    if (i >= n) return;
    const auto& x = inputs[i];
    if (!valid_input(x, options) || options.temperature_scale_K <= 0.0 ||
        options.residual_tolerance <= 0.0 || options.difference_step <= 0.0) {
        outputs[i] = failure(ImpellerBatchStatus::InvalidInput,0); return;
    }
    if (!(x.diffuser_outlet_radius_m>x.exit_diameter_m/2.0 &&
          x.diffuser_radial_steps>0 && x.exit_blockage==0.0)) {
        outputs[i] = failure(ImpellerBatchStatus::InvalidInput,0); return;
    }
    const auto inlet = thermomodel::detail::pr_mixture_evaluate_Tp_cuda(
        mixture,x.inlet_temperature_K,x.inlet_pressure_Pa);
    if (inlet.status != thermomodel::detail::PRMixtureCudaStatus::Success || inlet.rho <= 0.0) {
        outputs[i] = failure(ImpellerBatchStatus::InletThermoFailed,0); return;
    }

    double y[2]{x.initial_temperature_K/options.temperature_scale_K,x.initial_log_pressure};
    Residual current = evaluate_scaled(x,options,mixture,viscosity_model,inlet,y[0],y[1]);
    if (!current.valid) { outputs[i] = failure(ImpellerBatchStatus::OutletThermoFailed,0); return; }
    for (int iteration=0; iteration<=options.maximum_iterations; ++iteration) {
        if (current.norm <= options.residual_tolerance) {
            const DiffuserResult diffuser=integrate_diffuser(
                x,mixture,viscosity_model,options.use_iso20765_viscosity,current);
            if (!diffuser.valid) {
                outputs[i]=failure(ImpellerBatchStatus::OutletThermoFailed,iteration);
                return;
            }
            ImpellerBatchOutput out{};
            out.status=ImpellerBatchStatus::Converged; out.iterations=iteration;
            out.outlet_temperature_K=y[0]*options.temperature_scale_K;
            out.outlet_pressure_Pa=current.p; out.outlet_density_kg_m3=current.rho;
            out.meridional_velocity_m_s=current.vm; out.whirl_velocity_m_s=current.whirl;
            out.euler_work_J_kg=current.work; out.shaft_torque_N_m=current.torque;
            out.energy_residual_J_kg=current.er; out.entropy_residual_J_kgK=current.sr;
            out.residual_norm=current.norm; out.darcy_friction_factor=current.friction_factor;
            out.passage_entropy_generation_J_kgK=current.friction_entropy;
            out.diffuser_outlet_temperature_K=diffuser.T;
            out.diffuser_outlet_pressure_Pa=diffuser.p;
            out.diffuser_outlet_density_kg_m3=diffuser.rho;
            out.diffuser_outlet_meridional_velocity_m_s=diffuser.vm;
            out.diffuser_outlet_whirl_m_s=diffuser.whirl;
            out.diffuser_stagnation_enthalpy_residual_J_kg=diffuser.energy_residual;
            out.diffuser_entropy_change_J_kgK=diffuser.entropy_change;
            out.diffuser_integrated_entropy_J_kgK=diffuser.entropy_integral;
            out.diffuser_entropy_balance_residual_J_kgK=diffuser.entropy_residual;
            out.diffuser_angular_momentum_residual_m2_s=diffuser.angular_residual;
            out.diffuser_integrated_wall_torque_m2_s_per_kg=diffuser.wall_torque_integral;
            out.diffuser_mass_residual_kg_s=diffuser.mass_residual;
            outputs[i]=out; return;
        }
        if (iteration == options.maximum_iterations) break;

        double J[2][2]{};
        for (int j=0;j<2;++j) {
            const double h=options.difference_step*fmax(1.0,fabs(y[j]));
            double yp[2]{y[0],y[1]}, ym[2]{y[0],y[1]}; yp[j]+=h; ym[j]-=h;
            const Residual plus=evaluate_scaled(x,options,mixture,viscosity_model,inlet,yp[0],yp[1]);
            const Residual minus=evaluate_scaled(x,options,mixture,viscosity_model,inlet,ym[0],ym[1]);
            if (!plus.valid || !minus.valid) { outputs[i]=failure(ImpellerBatchStatus::JacobianFailed,iteration); return; }
            J[0][j]=(plus.f0-minus.f0)/(2.0*h);
            J[1][j]=(plus.f1-minus.f1)/(2.0*h);
        }
        const double n0=hypot(J[0][0],J[1][0]), n1=hypot(J[0][1],J[1][1]);
        if (!(n0>0.0 && n1>0.0)) { outputs[i]=failure(ImpellerBatchStatus::SingularJacobian,iteration); return; }
        const double a=J[0][0]/n0,c=J[1][0]/n0,b=J[0][1]/n1,d=J[1][1]/n1;
        const double det=a*d-b*c;
        if (!isfinite(det) || fabs(det)<=options.jacobian_singularity_tolerance) {
            outputs[i]=failure(ImpellerBatchStatus::SingularJacobian,iteration); return;
        }
        const double step0=(d*(-current.f0)-b*(-current.f1))/det/n0;
        const double step1=(-c*(-current.f0)+a*(-current.f1))/det/n1;
        double alpha=1.0; bool accepted=false;
        for (int backtrack=0;backtrack<=options.maximum_backtracks;++backtrack) {
            if (alpha<options.minimum_step) break;
            const Residual candidate=evaluate_scaled(x,options,mixture,viscosity_model,inlet,
                y[0]+alpha*step0,y[1]+alpha*step1);
            if (candidate.valid &&
                candidate.norm <= current.norm*(1.0-1.0e-4*alpha)) {
                y[0]+=alpha*step0; y[1]+=alpha*step1; current=candidate; accepted=true; break;
            }
            alpha*=0.5;
        }
        if (!accepted) { outputs[i]=failure(ImpellerBatchStatus::LineSearchFailed,iteration); return; }
    }
    outputs[i]=failure(ImpellerBatchStatus::MaximumIterations,options.maximum_iterations);
}

void check(cudaError_t e,const char* what)
{ if(e!=cudaSuccess) throw std::runtime_error(std::string(what)+": "+cudaGetErrorString(e)); }

template<class T> T* device_alloc(std::size_t n)
{ T* p=nullptr; check(cudaMalloc(reinterpret_cast<void**>(&p),n*sizeof(T)),"cudaMalloc"); return p; }

} // namespace

struct ImpellerBatchSolver::Impl {
    thermogpu::thermo::CUDAPRMixtureWorkspace thermo;
    transport::Iso20765Parameters viscosity_model;
    std::size_t capacity;
    ImpellerBatchInput* device_inputs = nullptr;
    ImpellerBatchOutput* device_outputs = nullptr;
    std::mutex solve_mutex;

    Impl(const thermomodel::PengRobinsonMixture& model,std::size_t n)
        : thermo(model,n),viscosity_model(transport::make_iso20765_parameters(model)),capacity(n)
    {
        device_inputs=device_alloc<ImpellerBatchInput>(capacity);
        try {
            device_outputs=device_alloc<ImpellerBatchOutput>(capacity);
        } catch (...) {
            cudaFree(device_inputs);
            device_inputs=nullptr;
            throw;
        }
    }

    ~Impl()
    {
        if (device_outputs) cudaFree(device_outputs);
        if (device_inputs) cudaFree(device_inputs);
    }
};

ImpellerBatchSolver::ImpellerBatchSolver(const thermomodel::PengRobinsonMixture& m,std::size_t n)
 : impl_(nullptr)
{
    if(n==0) throw std::invalid_argument("GPU impeller batch capacity must be positive");
    impl_=std::make_unique<Impl>(m,n);
}
ImpellerBatchSolver::~ImpellerBatchSolver()=default;
ImpellerBatchSolver::ImpellerBatchSolver(ImpellerBatchSolver&&) noexcept=default;
ImpellerBatchSolver& ImpellerBatchSolver::operator=(ImpellerBatchSolver&&) noexcept=default;
std::size_t ImpellerBatchSolver::capacity() const noexcept { return impl_?impl_->capacity:0; }

std::vector<ImpellerBatchOutput> ImpellerBatchSolver::solve(
    const std::vector<ImpellerBatchInput>& inputs,const ImpellerBatchOptions& options) const
{
    if(!impl_ || inputs.empty() || inputs.size()>impl_->capacity)
        throw std::invalid_argument("GPU impeller batch is empty or exceeds capacity");
    std::lock_guard<std::mutex> lock(impl_->solve_mutex);
    check(cudaMemcpy(impl_->device_inputs,inputs.data(),
        inputs.size()*sizeof(ImpellerBatchInput),cudaMemcpyHostToDevice),"copy inputs");
    constexpr unsigned threads=128;
    const unsigned blocks=static_cast<unsigned>((inputs.size()+threads-1)/threads);
    solve_kernel<<<blocks,threads>>>(impl_->device_inputs,impl_->device_outputs,
        inputs.size(),impl_->thermo.device_parameters(),impl_->viscosity_model,options);
    check(cudaGetLastError(),"launch impeller kernel");
    std::vector<ImpellerBatchOutput> outputs(inputs.size());
    // The blocking D2H copy waits for the default-stream kernel and surfaces
    // asynchronous execution errors, so a separate device-wide sync is redundant.
    check(cudaMemcpy(outputs.data(),impl_->device_outputs,
        outputs.size()*sizeof(ImpellerBatchOutput),cudaMemcpyDeviceToHost),"copy outputs");
    return outputs;
}

} // namespace turbocompgpu::cuda
