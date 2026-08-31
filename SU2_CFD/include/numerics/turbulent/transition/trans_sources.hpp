/*!
 * \file trans_sources.hpp
 * \brief Numerics classes for integration of source terms in transition problems.
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

#pragma once
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../scalar/scalar_sources.hpp"
#include "./trans_correlations.hpp"

/*!
 * \class CSourcePieceWise_TranLM
 * \brief Class for integrating the source terms of the LM transition model equations.
 * \ingroup SourceDiscr
 * \author S. Kang.
 */
template <class FlowIndices>
class CSourcePieceWise_TransLM final : public CNumerics {
 private:
  const FlowIndices idx; /*!< \brief Object to manage the access to the flow primitives. */

  const LM_ParsedOptions options;
  const bool axisymmetric;

  /*--- LM Closure constants ---*/
  const su2double c_e1 = 1.0;
  const su2double c_a1 = 2.0;
  const su2double c_e2 = 50.0;
  const su2double c_a2 = 0.06;
  const su2double sigmaf = 1.0;
  const su2double s1 = 2.0;
  const su2double c_theta = 0.03;
  const su2double c_CF = 0.6;
  const su2double sigmat = 2.0;

  TURB_FAMILY TurbFamily;
  su2double hRoughness;

  su2double IntermittencySep = 1.0;
  su2double IntermittencyEff = 1.0;

  su2double Residual[2];
  su2double* Jacobian_i[2];
  su2double Jacobian_Buffer[4];  // Static storage for the Jacobian (which needs to be pointer for return type).

  TransLMCorrelations TransCorrelations;

  /*!
   * \brief Add contribution from convection and diffusion due to axisymmetric formulation to 2D residual.
   */
  void ResidualAxisymmetricConvectionDiffusion() {
    if (Coord_i[1] < EPS) return;

    const su2double yinv = 1.0 / Coord_i[1];
    const su2double rhov = Density_i * V_i[idx.Velocity() + 1];
    const su2double diffusivity = Laminar_Viscosity_i + Eddy_Viscosity_i;

    const su2double cdgamma_axi = rhov * TransVar_i[0] - diffusivity * TransVar_Grad_i[0][1];
    const su2double cdretheta_axi = rhov * TransVar_i[1] - 2.0 * diffusivity * TransVar_Grad_i[1][1];

    Residual[0] -= yinv * Volume * cdgamma_axi;
    Residual[1] -= yinv * Volume * cdretheta_axi;

    const su2double velocity_radial = V_i[idx.Velocity() + 1];
    Jacobian_i[0][0] -= yinv * Volume * velocity_radial;
    Jacobian_i[0][1] -= 0.0;
    Jacobian_i[1][0] -= 0.0;
    Jacobian_i[1][1] -= yinv * Volume * velocity_radial;
  }

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_TransLM(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config)
      : CNumerics(val_nDim, 2, config), idx(val_nDim, config->GetnSpecies()), options(config->GetLMParsedOptions()), axisymmetric(config->GetAxisymmetric()){
    /*--- "Allocate" the Jacobian using the static buffer. ---*/
    Jacobian_i[0] = Jacobian_Buffer;
    Jacobian_i[1] = Jacobian_Buffer + 2;

    TurbFamily = TurbModelFamily(config->GetKind_Turb_Model());

    hRoughness = config->GethRoughness();

    TransCorrelations.SetOptions(options);

  }

  /*!
   * \brief Residual for source term integration.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override {
    /*--- ScalarVar[0] = k, ScalarVar[0] = w, TransVar[0] = gamma, and TransVar[0] = ReThetaT ---*/
    /*--- dU/dx = PrimVar_Grad[1][0] ---*/
    AD::StartPreacc();
    AD::SetPreaccIn(StrainMag_i);
    AD::SetPreaccIn(ScalarVar_i, nVar);
    AD::SetPreaccIn(ScalarVar_Grad_i, nVar, nDim);
    AD::SetPreaccIn(TransVar_i, nVar);
    AD::SetPreaccIn(TransVar_Grad_i, nVar, nDim);
    AD::SetPreaccIn(Volume);
    AD::SetPreaccIn(dist_i);
    AD::SetPreaccIn(&V_i[idx.Velocity()], nDim);
    AD::SetPreaccIn(PrimVar_Grad_i, nDim + idx.Velocity(), nDim);
    AD::SetPreaccIn(Vorticity_i, 3);
    if (axisymmetric) AD::SetPreaccIn(Coord_i[1]);

    su2double VorticityMag =
        sqrt(Vorticity_i[0] * Vorticity_i[0] + Vorticity_i[1] * Vorticity_i[1] + Vorticity_i[2] * Vorticity_i[2]);

    const su2double vel_u = V_i[idx.Velocity()];
    const su2double vel_v = V_i[1 + idx.Velocity()];
    const su2double vel_w = (nDim == 3) ? V_i[2 + idx.Velocity()] : 0.0;

    const su2double Velocity_Mag = sqrt(vel_u * vel_u + vel_v * vel_v + vel_w * vel_w);

    AD::SetPreaccIn(V_i[idx.Density()], V_i[idx.LaminarViscosity()], V_i[idx.EddyViscosity()]);

    Density_i = V_i[idx.Density()];
    Laminar_Viscosity_i = V_i[idx.LaminarViscosity()];
    Eddy_Viscosity_i = V_i[idx.EddyViscosity()];

    Residual[0] = 0.0;
    Residual[1] = 0.0;
    Jacobian_i[0][0] = 0.0;
    Jacobian_i[0][1] = 0.0;
    Jacobian_i[1][0] = 0.0;
    Jacobian_i[1][1] = 0.0;

    if (dist_i > 1e-10) {
      su2double Tu = 1.0;
      if (TurbFamily == TURB_FAMILY::KW) Tu = max(100.0 * sqrt(2.0 * ScalarVar_i[0] / 3.0) / Velocity_Mag, 0.027);
      if (TurbFamily == TURB_FAMILY::SA) Tu = config->GetTurbulenceIntensity_FreeStream() * 100;

      /*--- Corr_RetC correlation*/
      const su2double Corr_Rec = TransCorrelations.ReThetaC_Correlations(Tu, TransVar_i[1]);

      /*--- F_length correlation*/
      const su2double Corr_F_length = TransCorrelations.FLength_Correlations(Tu, TransVar_i[1]);

      /*--- F_length ---*/
      su2double F_length = 0.0;
      if (TurbFamily == TURB_FAMILY::KW) {
        const su2double r_omega = Density_i * dist_i * dist_i * ScalarVar_i[1] / Laminar_Viscosity_i;
        const su2double f_sub = exp(-pow(r_omega / 200.0, 2));
        F_length = Corr_F_length * (1. - f_sub) + 40.0 * f_sub;
        if (options.LMFAN) {
          F_length = 20.0 * (1.0 - f_sub) + 40.0 * f_sub;
        }
      }
      if (TurbFamily == TURB_FAMILY::SA) F_length = Corr_F_length;

      /*--- F_onset ---*/
      su2double R_t = 1.0;
      if (TurbFamily == TURB_FAMILY::KW) R_t = Density_i * ScalarVar_i[0] / Laminar_Viscosity_i / ScalarVar_i[1];
      if (TurbFamily == TURB_FAMILY::SA) R_t = Eddy_Viscosity_i / Laminar_Viscosity_i;

      const su2double Re_v = Density_i * dist_i * dist_i * StrainMag_i / Laminar_Viscosity_i;
      su2double F_onset1 = Re_v / (2.193 * Corr_Rec);
      su2double F_onset2 = 1.0;
      su2double F_onset3 = 1.0;

      /*--- Hypersonic Correction for streamwise & crossflow transition by Fan (Used in LMFAN) ---*/
      const su2double sos = V_i[idx.SoundSpeed()];
      const su2double rho_inf = config->GetDensity_FreeStream();
      const su2double p_inf = config->GetPressure_FreeStream();
      const su2double velU_inf = config->GetVelocity_FreeStream()[0];
      const su2double velV_inf = config->GetVelocity_FreeStream()[1];
      const su2double velW_inf = (nDim == 3) ? config->GetVelocity_FreeStream()[2] : 0.0;
      const su2double velMag_inf = pow(velU_inf * velU_inf + velV_inf * velV_inf + velW_inf * velW_inf, 0.5);
      const su2double gamma_Spec = config->GetGamma();
      const su2double p = V_i[idx.Pressure()];
      const su2double sos_inf = pow(config->GetTemperature_FreeStream() * config->GetGas_Constant() * gamma_Spec, 0.5);
      const su2double temperature_local = V_i[idx.Temperature()];
      const su2double Twall = config->GetFAN_Twall();

      su2double F_onset_s = 0.0;
      su2double F_onset_cf = 0.0;

      su2double rho_eL = 0.0, U_eL = 0.0, a_eL = 0.0, T_eL = 0.0, Ma_eL = 0.0, He = 0.0;

      rho_eL = pow(rho_inf, gamma_Spec) * p / p_inf;
      rho_eL = pow(rho_eL, 1 / gamma_Spec);
      U_eL = gamma_Spec / (gamma_Spec -1.0) * p_inf / rho_inf + 0.5 * velMag_inf * velMag_inf;
      U_eL -= gamma_Spec / (gamma_Spec -1.0) * p / rho_eL;
      U_eL = pow(U_eL * 2.0, 0.5);
      a_eL = sos_inf * sos_inf / (gamma_Spec -1) + velMag_inf * velMag_inf / 2.0;
      a_eL -= U_eL * U_eL / 2.0;
      a_eL = pow(a_eL * (gamma_Spec -1), 0.5);
      T_eL = a_eL * a_eL / gamma_Spec / config->GetGas_Constant();
      Ma_eL = U_eL / a_eL;

      su2double a1 = + 1.882e-4 * Ma_eL * Ma_eL * Ma_eL + 4.544e-3 * Ma_eL * Ma_eL - 1.954e-1 * Ma_eL + 1.748;
      su2double a2 = + 1.667e-4 * Ma_eL * Ma_eL * Ma_eL - 2.171e-3 * Ma_eL * Ma_eL - 2.937e-2 * Ma_eL - 0.5902;
      su2double a3 = - 8.928e-4 * Ma_eL * Ma_eL * Ma_eL + 2.041e-2 * Ma_eL * Ma_eL + 9.166e-2 * Ma_eL + 0.4975;

      su2double F_ratio = a1 * pow(T_eL / Twall, a2) + a3;

      if (nDim == 2) {
        He = 0.0;
      }
      else {
        su2double VelocityNormalized[3];
        VelocityNormalized[0] = vel_u / Velocity_Mag;
        VelocityNormalized[1] = vel_v / Velocity_Mag;
        VelocityNormalized[2] = vel_w / Velocity_Mag;

        su2double StreamwiseVort = 0.0;
        for (auto iDim =0u; iDim < nDim; iDim++) {
          StreamwiseVort += VelocityNormalized[iDim] * Vorticity_i[iDim];
        }
        StreamwiseVort = abs(StreamwiseVort);
        He = StreamwiseVort;
      }
      
      su2double delH_cf = 0.0;
      su2double H_cf = 0.0;
      su2double C_cf = 28.0;

      const su2double H_CF = He * dist_i / Velocity_Mag;
      const su2double delH_CF = H_CF *(1.0 + min(Eddy_Viscosity_i / Laminar_Viscosity_i, 0.4));
      H_cf = He * dist_i / Velocity_Mag;
      delH_cf = H_CF *(1.0 + min(Eddy_Viscosity_i / Laminar_Viscosity_i, 0.4));


      su2double F_Tu = config->GetTurbulenceIntensity_FreeStream();
      if (F_Tu < 0.001) C_cf = 45.0;

      if (options.LMFAN) {
        F_onset_s = Re_v / F_ratio / Corr_Rec;
        F_onset_cf = delH_cf * Re_v / F_ratio / C_cf;
        F_onset1 = max(F_onset_s, F_onset_cf);
      }

      if (TurbFamily == TURB_FAMILY::KW) {
        F_onset2 = min(max(F_onset1, pow(F_onset1, 4.0)), 2.0);
        F_onset3 = max(1.0 - pow(R_t / 2.5, 3.0), 0.0);
      }
      if (TurbFamily == TURB_FAMILY::SA) {
        F_onset2 = min(max(F_onset1, pow(F_onset1, 4.0)), 4.0);
        F_onset3 = max(2.0 - pow(R_t / 2.5, 3.0), 0.0);
      }
      const su2double F_onset = max(F_onset2 - F_onset3, 0.0);

      /*-- Gradient of velocity magnitude ---*/

      su2double dU_dx = 0.5 / Velocity_Mag * (2. * vel_u * PrimVar_Grad_i[1][0] + 2. * vel_v * PrimVar_Grad_i[2][0]);
      if (nDim == 3) dU_dx += 0.5 / Velocity_Mag * (2. * vel_w * PrimVar_Grad_i[3][0]);

      su2double dU_dy = 0.5 / Velocity_Mag * (2. * vel_u * PrimVar_Grad_i[1][1] + 2. * vel_v * PrimVar_Grad_i[2][1]);
      if (nDim == 3) dU_dy += 0.5 / Velocity_Mag * (2. * vel_w * PrimVar_Grad_i[3][1]);

      su2double dU_dz = 0.0;
      if (nDim == 3)
        dU_dz =
            0.5 / Velocity_Mag *
            (2. * vel_u * PrimVar_Grad_i[1][2] + 2. * vel_v * PrimVar_Grad_i[2][2] + 2. * vel_w * PrimVar_Grad_i[3][2]);

      su2double du_ds = vel_u / Velocity_Mag * dU_dx + vel_v / Velocity_Mag * dU_dy;
      if (nDim == 3) du_ds += vel_w / Velocity_Mag * dU_dz;

      /*-- Calculate blending function f_theta --*/
      su2double time_scale = 500.0 * Laminar_Viscosity_i / Density_i / Velocity_Mag / Velocity_Mag;
      if (options.LM2015)
        time_scale = min(time_scale,
                         Density_i * LocalGridLength_i * LocalGridLength_i / (Laminar_Viscosity_i + Eddy_Viscosity_i));
      const su2double theta_bl = TransVar_i[1] * Laminar_Viscosity_i / Density_i / Velocity_Mag;
      const su2double delta_bl = 7.5 * theta_bl;
      const su2double delta = 50.0 * VorticityMag * dist_i / Velocity_Mag * delta_bl + 1e-20;

      su2double f_wake = 0.0;
      if (TurbFamily == TURB_FAMILY::KW) {
        const su2double re_omega = Density_i * ScalarVar_i[1] * dist_i * dist_i / Laminar_Viscosity_i;
        f_wake = exp(-pow(re_omega / (1.0e+05), 2));
      }
      if (TurbFamily == TURB_FAMILY::SA) f_wake = 1.0;

      const su2double var1 = (TransVar_i[0] - 1.0 / c_e2) / (1.0 - 1.0 / c_e2);
      const su2double var2 = 1.0 - pow(var1, 2.0);
      const su2double f_theta = min(max(f_wake * exp(-pow(dist_i / delta, 4)), var2), 1.0);
      const su2double f_turb = exp(-pow(R_t / 4, 4));

      su2double f_theta_2 = 0.0;
      if (options.LM2015)
        f_theta_2 = min(f_wake * exp(-pow(dist_i / delta, 4.0)), 1.0);

      /*--- Corr_Ret correlation*/
      const su2double Corr_Ret_lim = 20.0;
      su2double f_lambda = 1.0;

      su2double Retheta_Error = 200.0, Retheta_old = 0.1;
      su2double lambda = 0.0;
      su2double Corr_Ret = 20.0;

      for (int iter = 0; iter < 100; iter++) {
        su2double theta = Corr_Ret * Laminar_Viscosity_i / Density_i / Velocity_Mag;
        lambda = Density_i * theta * theta / Laminar_Viscosity_i * du_ds;
        
        if (options.LMFAN) {
          lambda = lambda * (1.0 + (gamma_Spec - 1.0) / 2.0 * Ma_eL * Ma_eL);
        }

        lambda = min(max(-0.1, lambda), 0.1);

        if (lambda <= 0.0) {
          f_lambda = 1. - (-12.986 * lambda - 123.66 * lambda * lambda - 405.689 * lambda * lambda * lambda) *
                              exp(-pow(Tu / 1.5, 1.5));
        } else {
          f_lambda = 1. + 0.275 * (1. - exp(-35. * lambda)) * exp(-Tu / 0.5);
        }

        if (Tu <= 1.3) {
          Corr_Ret = f_lambda * (1173.51 - 589.428 * Tu + 0.2196 / Tu / Tu);
        } else {
          Corr_Ret = 331.5 * f_lambda * pow(Tu - 0.5658, -0.671);
        }

        if (options.LMFAN) {
          su2double fMa_eL = -83.16 * pow(Ma_eL, -4.095) + 1.509;
          fMa_eL = max(fMa_eL, 0.1);
          Corr_Ret = Corr_Ret / fMa_eL;
        }

        Corr_Ret = max(Corr_Ret, Corr_Ret_lim);

        Retheta_Error = fabs(Retheta_old - Corr_Ret) / Retheta_old;

        if (Retheta_Error < 0.0000001) {
          break;
        }

        Retheta_old = Corr_Ret;
      }

      /*-- Corr_RetT_SCF Correlations--*/
      su2double ReThetat_SCF = 0.0;
      if (options.LM2015) {
        su2double VelocityNormalized[3];
        VelocityNormalized[0] = vel_u / Velocity_Mag;
        VelocityNormalized[1] = vel_v / Velocity_Mag;
        if (nDim == 3) VelocityNormalized[2] = vel_w / Velocity_Mag;

        su2double StreamwiseVort = 0.0;
        for (auto iDim = 0u; iDim < nDim; iDim++) {
          StreamwiseVort += VelocityNormalized[iDim] * Vorticity_i[iDim];
        }
        StreamwiseVort = abs(StreamwiseVort);

        const su2double H_CF = StreamwiseVort * dist_i / Velocity_Mag;
        const su2double DeltaH_CF = H_CF * (1.0 + min(Eddy_Viscosity_i / Laminar_Viscosity_i, 0.4));
        const su2double DeltaH_CF_Minus = max(-1.0 * (0.1066 - DeltaH_CF), 0.0);
        const su2double DeltaH_CF_Plus = max(0.1066 - DeltaH_CF, 0.0);
        const su2double fDeltaH_CF_Minus = 75.0 * tanh(DeltaH_CF_Minus / 0.0125);
        const su2double fDeltaH_CF_Plus = 6200 * DeltaH_CF_Plus + 50000 * DeltaH_CF_Plus * DeltaH_CF_Plus;

        const su2double toll = 1e-5;
        su2double error = toll + 1.0;
        su2double thetat_SCF = 0.0;
        su2double rethetat_SCF_old = 20.0;
        const int nMax = 100;

        int iter;
        for (iter = 0; iter < nMax && error > toll; iter++) {
          thetat_SCF = rethetat_SCF_old * Laminar_Viscosity_i / (Density_i * (Velocity_Mag / 0.82));
          thetat_SCF = max(1e-20, thetat_SCF);

          ReThetat_SCF = -35.088 * log(hRoughness / thetat_SCF) + 319.51 + fDeltaH_CF_Plus - fDeltaH_CF_Minus;

          error = abs(ReThetat_SCF - rethetat_SCF_old) / rethetat_SCF_old;

          rethetat_SCF_old = ReThetat_SCF;
        }
      }

      /*-- production term of Intermeittency(Gamma) --*/
      const su2double Pg =
          F_length * c_a1 * Density_i * StrainMag_i * sqrt(F_onset * TransVar_i[0]) * (1.0 - c_e1 * TransVar_i[0]);

      /*-- destruction term of Intermeittency(Gamma) --*/
      const su2double Dg = c_a2 * Density_i * VorticityMag * TransVar_i[0] * f_turb * (c_e2 * TransVar_i[0] - 1.0);

      /*-- production term of ReThetaT --*/
      const su2double PRethetat = c_theta * Density_i / time_scale * (Corr_Ret - TransVar_i[1]) * (1.0 - f_theta);

      /*-- destruction term of ReThetaT --*/
      // It should not be with the minus sign but I put for consistency
      su2double DRethetat = 0.0;
      if (options.LM2015)
        DRethetat = -c_theta * (Density_i / time_scale) * c_CF * min(ReThetat_SCF - TransVar_i[1], 0.0) * f_theta_2;

      /*--- Source ---*/
      Residual[0] += (Pg - Dg) * Volume;
      Residual[1] += (PRethetat - DRethetat) * Volume;

      /*--- Implicit part ---*/
      Jacobian_i[0][0] = (F_length * c_a1 * StrainMag_i * sqrt(F_onset) *
                              (0.5 * pow(TransVar_i[0], -0.5) - 1.5 * c_e1 * pow(TransVar_i[0], 0.5)) -
                          c_a2 * VorticityMag * f_turb * (2.0 * c_e2 * TransVar_i[0] - 1.0)) *
                         Volume;
      Jacobian_i[0][1] = 0.0;
      Jacobian_i[1][0] = 0.0;
      Jacobian_i[1][1] = -c_theta / time_scale * (1.0 - f_theta) * Volume;
      if (options.LM2015 && ReThetat_SCF - TransVar_i[1] < 0)
        Jacobian_i[1][1] += (c_theta / time_scale) * c_CF * f_theta_2 * Volume;
    }

    if (axisymmetric) ResidualAxisymmetricConvectionDiffusion();

    AD::SetPreaccOut(Residual, nVar);
    AD::EndPreacc();

    return ResidualType<>(Residual, Jacobian_i, nullptr);
  }
};
