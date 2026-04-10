/*!
 * \file CTransLMSolver.cpp
 * \brief Main subroutines for Langtry-Menter Transition model solver.
 * \author A. Aranake, S. Kang.
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../include/solvers/CTransLMSolver.hpp"
#include "../../include/variables/CTransLMVariable.hpp"
#include "../../include/variables/CFlowVariable.hpp"
#include "../../include/variables/CTurbSAVariable.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

/*---  This is the implementation of the Langtry-Menter transition model.
       The main reference for this model is:Langtry, Menter, AIAA J. 47(12) 2009
       DOI: https://doi.org/10.2514/1.42362 ---*/

// Note: TransLM seems to use rho*gamma, rho*Re_sigma as Solution variables, thus Conservative=true

CTransLMSolver::CTransLMSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh)
    : CTurbSolver(geometry, config, true) {
  unsigned long iPoint;
  ifstream restart_file;
  string text_line;

  bool multizone = config->GetMultizone_Problem();

  /*--- Dimension of the problem --> 2 Transport equations (intermittency, Reth) ---*/
  nVar = 2;
  nPrimVar = 2;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();

  /*--- Define variables needed for transition from config file */
  options = config->GetLMParsedOptions();
  TransCorrelations.SetOptions(options);
  TurbFamily = TurbModelFamily(config->GetKind_Turb_Model());

  /*--- Single grid simulation ---*/

  if (iMesh == MESH_0) {

    /*--- Define some auxiliary vector related with the residual ---*/

    Residual_RMS.resize(nVar,0.0);
    Residual_Max.resize(nVar,0.0);
    Point_Max.resize(nVar,0);
    Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

    /*--- Initialization of the structure of the whole Jacobian ---*/

    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (LM transition model)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config, ReducerStrategy);
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    System.SetxIsZero(true);

    if (ReducerStrategy)
      EdgeFluxes.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);

    /*--- Initialize the BGS residuals in multizone problems. ---*/
    if (multizone){
      Residual_BGS.resize(nVar,0.0);
      Residual_Max_BGS.resize(nVar,0.0);
      Point_Max_BGS.resize(nVar,0);
      Point_Max_Coord_BGS.resize(nVar,nDim) = su2double(0.0);
    }

  }

  /*--- Initialize lower and upper limits---*/
  lowerlimit[0] = 1.0e-4;
  upperlimit[0] = 5.0;

  lowerlimit[1] = 1.0e-4;
  upperlimit[1] = 1.0e15;

  /*--- Far-field flow state quantities and initialization. ---*/
  const su2double Intensity = config->GetTurbulenceIntensity_FreeStream()*100.0;

  const su2double Intermittency_Inf  = 1.0;
  su2double ReThetaT_Inf = 100.0;

  /*--- Momentum thickness Reynolds number, initialized from freestream turbulent intensity*/
  if (Intensity <= 1.3) {
    if(Intensity >=0.027) {
      ReThetaT_Inf = (1173.51-589.428*Intensity+0.2196/(Intensity*Intensity));
    }
    else {
      ReThetaT_Inf = (1173.51-589.428*Intensity+0.2196/(0.27*0.27));
    }
  }
  else if(Intensity>1.3) {
    ReThetaT_Inf = 331.5*pow(Intensity-0.5658,-0.671);
  }

  Solution_Inf[0] = Intermittency_Inf;
  Solution_Inf[1] = ReThetaT_Inf;

  /*--- Initialize the solution to the far-field state everywhere. ---*/
  nodes = new CTransLMVariable(Intermittency_Inf, ReThetaT_Inf, 1.0, 1.0, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, MPI_QUANTITIES::SOLUTION);
  CompleteComms(geometry, config, MPI_QUANTITIES::SOLUTION);

  /*--- Initializate quantities for SlidingMesh Interface ---*/

  SlidingState.resize(nMarker);
  SlidingStateNodes.resize(nMarker);

  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {
      SlidingState[iMarker].resize(nVertex[iMarker], nPrimVar+1) = nullptr;
      SlidingStateNodes[iMarker].resize(nVertex[iMarker],0);
    }
  }

  /*-- Allocation of inlets has to happen in derived classes (not CTurbSolver),
    due to arbitrary number of turbulence variables ---*/

  Inlet_TurbVars.resize(nMarker);
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_TurbVars[iMarker].resize(nVertex[iMarker],nVar);
    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; ++iVertex) {
      Inlet_TurbVars[iMarker](iVertex,0) = Intermittency_Inf;
      Inlet_TurbVars[iMarker](iVertex,1) = ReThetaT_Inf;
    }
  }

   const su2double CFL = config->GetCFL(MGLevel)*config->GetCFLRedCoeff_Turb();
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    nodes->SetLocalCFL(iPoint, CFL);
  }
  Min_CFL_Local = CFL;
  Max_CFL_Local = CFL;
  Avg_CFL_Local = CFL;

  /*--- Add the solver name. ---*/
  SolverName = "LM model";

}

void CTransLMSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
         unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  SU2_OMP_SAFE_GLOBAL_ACCESS(config->SetGlobalParam(config->GetKind_Solver(), RunTime_EqSystem);)

  /*--- Upwind second order reconstruction and gradients ---*/
  CommonPreprocessing(geometry, config, Output);
}

void CTransLMSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {

  /*--- Compute LM model gradients. ---*/

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetSolution_Gradient_GG(geometry, config, -1);
  }
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetSolution_Gradient_LS(geometry, config, -1);
  }

  AD::StartNoSharedReading();
  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());
  auto* turbNodes = su2staticcast_p<CTurbVariable*>(solver_container[TURB_SOL]->GetNodes());

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint ++) {

    // Here the nodes already have the new solution, thus I have to compute everything from scratch

    const su2double rho = flowNodes->GetDensity(iPoint);
    const su2double mu  = flowNodes->GetLaminarViscosity(iPoint);
    const su2double muT = turbNodes->GetmuT(iPoint);
    const su2double dist = geometry->nodes->GetWall_Distance(iPoint);
    su2double VorticityMag = GeometryToolbox::Norm(3, flowNodes->GetVorticity(iPoint));
    su2double StrainMag =flowNodes->GetStrainMag(iPoint);
    VorticityMag = max(VorticityMag, 1e-12);
    StrainMag = max(StrainMag, 1e-12); // safety against division by zero
    const su2double Intermittency = nodes->GetSolution(iPoint,0);
    const su2double Re_t = nodes->GetSolution(iPoint,1);
    const su2double Re_v = rho*dist*dist*StrainMag/mu;
    const su2double vel_u = flowNodes->GetVelocity(iPoint, 0);
    const su2double vel_v = flowNodes->GetVelocity(iPoint, 1);
    const su2double vel_w = (nDim ==3) ? flowNodes->GetVelocity(iPoint, 2) : 0.0;
    const su2double VelocityMag = sqrt(vel_u*vel_u + vel_v*vel_v + vel_w*vel_w);
    su2double omega = 0.0;
    su2double k = 0.0;
    if(TurbFamily == TURB_FAMILY::KW){
      omega = turbNodes->GetSolution(iPoint,1);
      k = turbNodes->GetSolution(iPoint,0);
    }
    su2double Tu = 1.0;
    if(TurbFamily == TURB_FAMILY::KW)
      Tu = max(100.0*sqrt( 2.0 * k / 3.0 ) / VelocityMag,0.027);
    if(TurbFamily == TURB_FAMILY::SA)
      Tu = config->GetTurbulenceIntensity_FreeStream()*100;

    const su2double Corr_Rec = TransCorrelations.ReThetaC_Correlations(Tu, Re_t);

    su2double R_t = 1.0;
    if(TurbFamily == TURB_FAMILY::KW)
      R_t = rho*k/ mu/ omega;
    if(TurbFamily == TURB_FAMILY::SA)
      R_t = muT/ mu;

    const su2double f_reattach = exp(-pow(R_t/20,4));
    su2double f_wake = 0.0;
    if(TurbFamily == TURB_FAMILY::KW){
      const su2double re_omega = rho*omega*dist*dist/mu;
      f_wake = exp(-pow(re_omega/(1.0e+05),2));
    }
    if(TurbFamily == TURB_FAMILY::SA)
      f_wake = 1.0;

    const su2double theta_bl   = Re_t*mu / rho /VelocityMag;
    const su2double delta_bl   = 7.5*theta_bl;
    const su2double delta      = 50.0*VorticityMag*dist/VelocityMag*delta_bl + 1e-20;
    const su2double var1 = (Intermittency-1.0/50.0)/(1.0-1.0/50.0);
    const su2double var2 = 1.0 - pow(var1,2.0);
    const su2double f_theta = min(max(f_wake*exp(-pow(dist/delta, 4)), var2), 1.0);
    su2double Intermittency_Sep = 2.0*max(0.0, Re_v/(3.235*Corr_Rec)-1.0)*f_reattach;
    Intermittency_Sep = min(Intermittency_Sep,2.0)*f_theta;
    Intermittency_Sep = min(max(0.0, Intermittency_Sep), 2.0);
    nodes -> SetIntermittencySep(iPoint, Intermittency_Sep);
    nodes -> SetIntermittencyEff(iPoint, Intermittency_Sep);

    if(config->GetLMParsedOptions().LMFAN) {
      const su2double sos = flowNodes->GetSoundSpeed(iPoint);
      const su2double rho_inf = config->GetDensity_FreeStream();
      const su2double p_inf = config->GetPressure_FreeStream();
      const su2double velU_inf = config->GetVelocity_FreeStream()[0];
      const su2double velV_inf = config->GetVelocity_FreeStream()[1];
      const su2double velW_inf = (nDim == 3) ? config->GetVelocity_FreeStream()[3] : 0.0;
      const su2double velMag_inf = pow(velU_inf * velU_inf + velV_inf * velV_inf + velW_inf * velW_inf, 0.5);
      const su2double gamma_Spec = config->GetGamma();
      const su2double p = flowNodes->GetPressure(iPoint);
      const su2double sos_inf = pow(config->GetTemperature_FreeStream() * config->GetGas_Constant() * gamma_Spec, 0.5);
      const su2double temperature_local = flowNodes->GetTemperature(iPoint);
      const su2double Twall = config->GetFAN_Twall();

      const su2double M_inf = velMag_inf / sos_inf;
      const su2double T0 = config->GetTemperature_FreeStream() * (1 + (gamma_Spec - 1.0) / 2.0 * M_inf * M_inf);
      const su2double dist_i = geometry->nodes->GetWall_Distance(iPoint);
      const su2double StrainMag_i = flowNodes->GetStrainMag(iPoint);
      const su2double Eddy_Viscosity_i = turbNodes->GetmuT(iPoint);
      const su2double Laminar_Viscosity_i = flowNodes->GetLaminarViscosity(iPoint);
      const su2double Density_i = flowNodes->GetDensity(iPoint);
      const su2double Volum_i = geometry->nodes->GetVolume(iPoint);
      const su2double cordix = geometry->nodes->GetCoord(iPoint, 0);
      const su2double cordiy = geometry->nodes->GetCoord(iPoint, 1);
      const su2double turb_k = turbNodes->GetSolution(iPoint, 0);
      const su2double turb_w = turbNodes->GetSolution(iPoint, 0);

      su2double rho_eL = 0.0;
      su2double U_eL = 0.0;
      su2double a_eL = 0.0;
      su2double T_eL = 0.0;
      su2double M_eL = 0.0;
      su2double He = 0.0;
      su2double VorticityMag = 0.0;

      rho_eL = pow(rho_inf, gamma_Spec) * p / p_inf;
      rho_eL = pow(rho_eL, 1 / gamma_Spec);
      U_eL = gamma_Spec / (gamma_Spec -1.0) * p_inf / rho_inf + 0.5 * velMag_inf * velMag_inf;
      U_eL -= gamma_Spec / (gamma_Spec -1.0) * p / rho_eL;
      U_eL = pow(U_eL * 2.0, 0.5);
      a_eL = sos_inf * sos_inf / (gamma_Spec -1) + velMag_inf * velMag_inf / 2.0;
      a_eL -= U_eL * U_eL / 2.0;
      a_eL = pow(a_eL * (gamma_Spec -1), 0.5);
      T_eL = a_eL * a_eL / gamma_Spec / config->GetGas_Constant();
      M_eL = U_eL / a_eL;

      const su2double mu_eL = 0.00001716 * pow(T_eL / 273.15, 1.5) * (273.15 + 110.4) / (T_eL + 110.4);

      su2double F_onset1 = Re_v / (2.193 * Corr_Rec);
      su2double F_onset2 = 1.0;
      su2double F_onset3 = 1.0;
      su2double F_onset_s = 0.0;
      su2double F_onset_cf = 0.0;
      su2double StreamwiseVort = 0.0;

      if(nDim = 2) {
        VorticityMag = sqrt(flowNodes->GetVorticity(iPoint)[0] * flowNodes->GetVorticity(iPoint)[0] + flowNodes->GetVorticity(iPoint)[1] * flowNodes->GetVorticity(iPoint)[1]);
      }
      else {
        su2double VelocityNormalized[3];
        VelocityNormalized[0] = vel_u / VelocityMag;
        VelocityNormalized[1] = vel_v / VelocityMag;
        VelocityNormalized[2] = vel_w / VelocityMag;
        for (auto iDim =0u; iDim < nDim; iDim++) {
          StreamwiseVort += VelocityNormalized[iDim] * flowNodes->GetVorticity(iPoint)[iDim];
        }
        StreamwiseVort = abs(StreamwiseVort);
        
        const su2double unitU =  vel_u / VelocityMag;
        const su2double unitV =  vel_v / VelocityMag;
        const su2double unitW =  vel_w / VelocityMag;

        const su2double vorticity_x = flowNodes->GetVorticity(iPoint)[0];
        const su2double vorticity_y = flowNodes->GetVorticity(iPoint)[1];
        const su2double vorticity_z = flowNodes->GetVorticity(iPoint)[2];

        const su2double UVor_x = unitU * vorticity_x;
        const su2double VVor_y = unitV * vorticity_y;
        const su2double WVor_z = unitW * vorticity_z;

        VorticityMag = (vorticity_x * vorticity_x + vorticity_y * vorticity_y + vorticity_z * vorticity_z);
        He = pow(UVor_x * UVor_x + VVor_y * VVor_y + WVor_z * WVor_z, 0.5);
      }

      su2double a1 = + 1.882e-4 * M_eL * M_eL * M_eL + 4.544e-3 * M_eL * M_eL - 1.954e-1 * M_eL + 1.748;
      su2double a2 = + 1.667e-4 * M_eL * M_eL * M_eL - 2.171e-3 * M_eL * M_eL - 2.937e-2 * M_eL - 0.5902;
      su2double a3 = - 8.928e-4 * M_eL * M_eL * M_eL + 2.041e-2 * M_eL * M_eL + 9.166e-2 * M_eL + 0.4975;

      su2double F_ratio = a1 * pow(T_eL / Twall, a2) + a3;

      su2double delH_cf = 0.0;
      su2double H_cf = 0.0;
      su2double C_cf = 28.0;

      const su2double H_CF = He * dist_i / VelocityMag;
      const su2double delH_CF = H_CF *(1.0 + min(Eddy_Viscosity_i / Laminar_Viscosity_i, 0.4));
      H_cf = He * dist_i / VelocityMag;
      delH_cf = H_CF *(1.0 + min(Eddy_Viscosity_i / Laminar_Viscosity_i, 0.4));

      su2double F_Tu = config->GetTke_FreeStream();
      F_Tu = config->GetTurbulenceIntensity_FreeStream();
      if (F_Tu < 0.001) C_cf = 45.0;

      F_onset_s = Re_v / F_ratio / Corr_Rec;
      F_onset_cf = delH_cf * Re_v / F_ratio / C_cf;
      F_onset1 = max(F_onset_s, F_onset_cf);
      F_onset2 = min(max(F_onset1, pow(F_onset1, 4.0)), 2.0);
      F_onset3 = max(1.0 - pow(R_t / 2.5, 3.0), 0.0);

      const su2double F_onset = max(F_onset2 - F_onset3, 0.0);

      if (dist_i < 1e-10) {
        su2double StreamwiseVort_temp = 0.0;
        su2double He_temp = 0.0;
        su2double H_CF_temp = 0.0;
        su2double delH_cf_temp = 0.0;
        su2double F_onset_cf_temp = 0.0;
    
        nodes->SetLM_Wonder_Func(iPoint, M_eL, T_eL, rho_eL, StreamwiseVort_temp, He_temp, H_CF_temp, delH_cf_temp, F_onset_s, F_onset_cf_temp, F_onset);
      }
      else {
        nodes->SetLM_Wonder_Func(iPoint, M_eL, T_eL, rho_eL, StreamwiseVort, He, H_CF, delH_cf, F_onset_s, F_onset_cf, F_onset);
      }
    }

  }
  END_SU2_OMP_FOR

  AD::EndNoSharedReading();
}


void CTransLMSolver::Viscous_Residual(const unsigned long iEdge, const CGeometry* geometry, CSolver** solver_container,
                                     CNumerics* numerics, const CConfig* config) {

  /*--- Define an object to set solver specific numerics contribution. ---*/

  auto SolverSpecificNumerics = [&](unsigned long iPoint, unsigned long jPoint) {};

  /*--- Now instantiate the generic implementation with the functor above. ---*/

  Viscous_Residual_impl(SolverSpecificNumerics, iEdge, geometry, solver_container, numerics, config);
}


void CTransLMSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                     CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());
  //auto* turbNodes = su2staticcast_p<CFlowVariable*>(solver_container[TURB_SOL]->GetNodes());
  CVariable* turbNodes = solver_container[TURB_SOL]->GetNodes();

  /*--- Pick one numerics object per thread. ---*/
  auto* numerics = numerics_container[SOURCE_FIRST_TERM + omp_get_thread_num()*MAX_TERMS];

  /*--- Loop over all points. ---*/

  AD::StartNoSharedReading();

  SU2_OMP_FOR_DYN(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {



    /*--- Conservative variables w/o reconstruction ---*/

    numerics->SetPrimitive(flowNodes->GetPrimitive(iPoint), nullptr);

    /*--- Gradient of the primitive and conservative variables ---*/

    numerics->SetPrimVarGradient(flowNodes->GetGradient_Primitive(iPoint), nullptr);

    /*--- Turbulent variables w/o reconstruction, and its gradient ---*/
    /*--- ScalarVar & ScalarVarGradient : Turbulence model solution(k&w) ---*/

    numerics->SetScalarVar(turbNodes->GetSolution(iPoint), nullptr);
    numerics->SetScalarVarGradient(turbNodes->GetGradient(iPoint), nullptr);

    /*--- Transition variables w/o reconstruction, and its gradient ---*/

    numerics->SetTransVar(nodes->GetSolution(iPoint), nullptr);
    numerics->SetTransVarGradient(nodes->GetGradient(iPoint), nullptr);

    /*--- Set volume ---*/

    numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

    /*--- Set distance to the surface ---*/

    numerics->SetDistance(geometry->nodes->GetWall_Distance(iPoint), 0.0);

    /*--- Set vorticity and strain rate magnitude ---*/

    numerics->SetVorticity(flowNodes->GetVorticity(iPoint), nullptr);

    numerics->SetStrainMag(flowNodes->GetStrainMag(iPoint), 0.0);

    /*--- Set coordinate (for debugging) ---*/
    numerics->SetCoord(geometry->nodes->GetCoord(iPoint), nullptr);

    if (options.LM2015) {
      /*--- Set local grid length (for LM2015)*/
      numerics->SetLocalGridLength(geometry->nodes->GetMaxLength(iPoint));
    }

    /*--- Compute the source term ---*/

    auto residual = numerics->ComputeResidual(config);

    /*--- Subtract residual and the Jacobian ---*/

    LinSysRes.SubtractBlock(iPoint, residual);
    if (implicit) Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

  }
  END_SU2_OMP_FOR

  AD::EndNoSharedReading();

}

void CTransLMSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                       CConfig *config, unsigned short iMesh) {
}

void CTransLMSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                      CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Allocate the value at the infinity ---*/

      auto V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      conv_numerics->SetPrimitive(V_domain, V_infty);

      /*--- Set turbulent variable at the wall, and at infinity ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), Solution_Inf);

      /*--- Set Normal (it is necessary to change the sign) ---*/
      /*--- It's mean wall normal zero flux. */

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      /*--- Grid Movement ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute residuals and Jacobians ---*/

      auto residual = conv_numerics->ComputeResidual(config);

      /*--- Add residuals and Jacobians ---*/

      LinSysRes.AddBlock(iPoint, residual);
      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
    }
  }
  END_SU2_OMP_FOR

}

void CTransLMSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                        CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  BC_HeatFlux_Wall(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);

}

void CTransLMSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                unsigned short val_marker) {
  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      /*--- Allocate the value at the inlet ---*/

      auto V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Non-dimensionalize Inlet_TurbVars if Inlet-Files are used. ---*/
      su2double Inlet_Vars[MAXNVAR];
      Inlet_Vars[0] = Inlet_TurbVars[val_marker][iVertex][0];
      Inlet_Vars[1] = Inlet_TurbVars[val_marker][iVertex][1];
      if (config->GetInlet_Profile_From_File()) {
        Inlet_Vars[0] /= pow(config->GetVelocity_Ref(), 2);
        Inlet_Vars[1] *= config->GetViscosity_Ref() / (config->GetDensity_Ref() * pow(config->GetVelocity_Ref(), 2));
      }

      /*--- Set the LM variable states. ---*/
      /*--- Load the inlet transition LM model variables (uniform by default). ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), Inlet_Vars);

      /*--- Set various other quantities in the solver class ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
    }
  }
  END_SU2_OMP_FOR

}

void CTransLMSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                               CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  BC_Far_Field(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
}

void CTransLMSolver::LoadRestart(CGeometry** geometry, CSolver*** solver, CConfig* config, int val_iter,
                                  bool val_update_geo) {

  const string restart_filename = config->GetFilename(config->GetSolution_FileName(), "", val_iter);

  /*--- To make this routine safe to call in parallel most of it can only be executed by one thread. ---*/
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

    if (config->GetRead_Binary_Restart()) {
      Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
    } else {
      Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
    }

    /*--- Skip flow variables ---*/

    unsigned short skipVars = nDim + solver[MESH_0][FLOW_SOL]->GetnVar() + solver[MESH_0][TURB_SOL] ->GetnVar();

    /*--- Adjust the number of solution variables in the incompressible
     restart. We always carry a space in nVar for the energy equation in the
     mean flow solver, but we only write it to the restart if it is active.
     Therefore, we must reduce skipVars here if energy is inactive so that
     the turbulent variables are read correctly. ---*/

    const bool incompressible = (config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE);
    const bool energy = config->GetEnergy_Equation();
    const bool weakly_coupled_heat = config->GetWeakly_Coupled_Heat();

    if (incompressible && ((!energy) && (!weakly_coupled_heat))) skipVars--;

    /*--- Load data from the restart into correct containers. ---*/

    unsigned long counter = 0;
    for (auto iPoint_Global = 0ul; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++) {
      /*--- Retrieve local index. If this node from the restart file lives
       on the current processor, we will load and instantiate the vars. ---*/

      const auto iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

      if (iPoint_Local > -1) {
        /*--- We need to store this point's data, so jump to the correct
         offset in the buffer of data from the restart file and load it. ---*/

        const auto index = counter * Restart_Vars[1] + skipVars;
        for (auto iVar = 0u; iVar < nVar; iVar++) nodes->SetSolution(iPoint_Local, iVar, Restart_Data[index + iVar]);
        nodes ->SetIntermittencySep(iPoint_Local,  Restart_Data[index + 2]);
        nodes ->SetIntermittencyEff(iPoint_Local,  Restart_Data[index + 3]);
        
        if(config->GetLMParsedOptions().LMFAN) {
          nodes ->SetLM_Wonder_Func(iPoint_Local, Restart_Data[index + 4], Restart_Data[index + 5], Restart_Data[index + 6],
                                    Restart_Data[index + 7], Restart_Data[index + 8], Restart_Data[index + 9],
                                    Restart_Data[index + 10], Restart_Data[index + 11], Restart_Data[index + 12],
                                    Restart_Data[index + 13]);
        }

        /*--- Increment the overall counter for how many points have been loaded. ---*/
        counter++;
      }
    }

    /*--- Detect a wrong solution file ---*/

    if (counter != nPointDomain) {
      SU2_MPI::Error(string("The solution file ") + restart_filename + string(" does not match with the mesh file!\n") +
                         string("This can be caused by empty lines at the end of the file."),
                     CURRENT_FUNCTION);
    }

  }  // end SU2_OMP_MASTER, pre and postprocessing are thread-safe.
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  /*--- MPI solution and compute the eddy viscosity ---*/

  solver[MESH_0][TRANS_SOL]->InitiateComms(geometry[MESH_0], config, MPI_QUANTITIES::SOLUTION);
  solver[MESH_0][TRANS_SOL]->CompleteComms(geometry[MESH_0], config, MPI_QUANTITIES::SOLUTION);

  /*--- For turbulent+species simulations the solver Pre-/Postprocessing is done by the species solver. ---*/
  if (config->GetKind_Species_Model() == SPECIES_MODEL::NONE) {
    solver[MESH_0][FLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER,
                                            RUNTIME_FLOW_SYS, false);
    solver[MESH_0][TURB_SOL]->Postprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0);
    solver[MESH_0][TRANS_SOL]->Postprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0);
  }

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/

  for (auto iMesh = 1u; iMesh <= config->GetnMGLevels(); iMesh++) {

    MultigridRestriction(*geometry[iMesh - 1], solver[iMesh - 1][TRANS_SOL]->GetNodes()->GetSolution(),
                         *geometry[iMesh], solver[iMesh][TRANS_SOL]->GetNodes()->GetSolution());
    solver[iMesh][TRANS_SOL]->InitiateComms(geometry[iMesh], config, MPI_QUANTITIES::SOLUTION);
    solver[iMesh][TRANS_SOL]->CompleteComms(geometry[iMesh], config, MPI_QUANTITIES::SOLUTION);

    if (config->GetKind_Species_Model() == SPECIES_MODEL::NONE) {
      solver[iMesh][FLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS,
                                            false);
      solver[iMesh][TRANS_SOL]->Postprocessing(geometry[iMesh], solver[iMesh], config, iMesh);
    }
  }

  /*--- Go back to single threaded execution. ---*/
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    /*--- Delete the class memory that is used to load the restart. ---*/

    Restart_Vars.clear();
    Restart_Data.clear();
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

}