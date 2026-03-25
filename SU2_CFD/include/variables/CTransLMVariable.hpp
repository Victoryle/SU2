/*!
 * \file CTransLMVariable.hpp
 * \brief Declaration of the variables of the transition model.
 * \author F. Palacios, T. Economon
 * \version 7.5.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "CTurbVariable.hpp"

/*!
 * \class CTransLMVariable
 * \brief Transition model variables.
 * \ingroup Turbulence_Model
 * \author A. Bueno, S. Kang.
 */

class CTransLMVariable final : public CTurbVariable {
protected:
  VectorType Intermittency_Eff;
  VectorType Intermittency_Sep;
  VectorType TempVar1, TempVar2, TempVar3, TempVar4, TempVar5, TempVar6;
  VectorType TempVar7, TempVar8, TempVar9, TempVar10; 

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] Intermittency - intermittency(gamma) (initialization value).
   * \param[in] ReThetaT - momentum thickness Reynolds number(ReThetaT)(initialization value).
   * \param[in] gammaSep - separation intermittency(gamma) (initialization value).
   * \param[in] gammaEff - effective intermittency(gamma) (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CTransLMVariable(su2double Intermittency, su2double ReThetaT, su2double gammaSep, su2double gammaEff, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTransLMVariable() override = default;

  /*!
   * \brief Set Separation intermittency.
   */
  void SetIntermittencySep(unsigned long iPoint, su2double val_Intermittency_sep) override;

  /*!
   * \brief Set Effective intermittency.
   */
  void SetIntermittencyEff(unsigned long iPoint, su2double val_Intermittency_sep) override;

  /*!
   * \brief Set Wonder variable.
   */
  void SetLM_Wonder_Func(unsigned long iPoint, su2double tempVar1, su2double tempVar2, su2double tempVar3,
                                      su2double tempVar4, su2double tempVar5, su2double tempVar6,su2double tempVar7,
                                      su2double tempVar8, su2double tempVar9, su2double tempVar10) override;

  /*!
   * \brief Value of Wonder variable.
   */
  inline su2double GetLM_Wonder_Func_var1(unsigned long iPoint) const override { return TempVar1(iPoint);}
  inline su2double GetLM_Wonder_Func_var2(unsigned long iPoint) const override { return TempVar2(iPoint);}
  inline su2double GetLM_Wonder_Func_var3(unsigned long iPoint) const override { return TempVar3(iPoint);}
  inline su2double GetLM_Wonder_Func_var4(unsigned long iPoint) const override { return TempVar4(iPoint);}
  inline su2double GetLM_Wonder_Func_var5(unsigned long iPoint) const override { return TempVar5(iPoint);}
  inline su2double GetLM_Wonder_Func_var6(unsigned long iPoint) const override { return TempVar6(iPoint);}
  inline su2double GetLM_Wonder_Func_var7(unsigned long iPoint) const override { return TempVar7(iPoint);}
  inline su2double GetLM_Wonder_Func_var8(unsigned long iPoint) const override { return TempVar8(iPoint);}
  inline su2double GetLM_Wonder_Func_var9(unsigned long iPoint) const override { return TempVar9(iPoint);}
  inline su2double GetLM_Wonder_Func_var10(unsigned long iPoint) const override { return TempVar10(iPoint);}

  /*!
   * \brief Calculate effective intermittency.
   */
  inline su2double GetIntermittencyEff(unsigned long iPoint) const override { return Intermittency_Eff(iPoint); }

  /*!
   * \brief Value of separation intermittency.
   */
  inline su2double GetIntermittencySep(unsigned long iPoint) const override { return Intermittency_Sep(iPoint); }

};
