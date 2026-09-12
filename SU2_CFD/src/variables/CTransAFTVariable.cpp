/*!
 * \file CTransLMVariable.cpp
 * \brief Definition of the solution fields.
 * \author A. Aranake, S. Kang
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


#include "../../include/variables/CTransAFTVariable.hpp"

CTransAFTVariable::CTransAFTVariable(su2double AF1, su2double AF2, su2double gammaAlgeEff, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config)
  : CTurbVariable(npoint, ndim, nvar, config) {

  for(unsigned long iPoint=0; iPoint<nPoint; ++iPoint)
  {
    Solution(iPoint,0) = AF1;
    Solution(iPoint,1) = AF2;
  }

  Solution_Old = Solution;

  /*--- Setting CTransLMVariable of intermittency_Eff---*/
  Intermittency_Alge_Eff.resize(nPoint) = gammaAlgeEff;

}

void CTransAFTVariable::SetIntermittencyAlgeEff(unsigned long iPoint, su2double val_Intermittency_Alge_Eff) {

  /*--- Effective intermittency ---*/
  Intermittency_Alge_Eff(iPoint) = val_Intermittency_Alge_Eff;

}