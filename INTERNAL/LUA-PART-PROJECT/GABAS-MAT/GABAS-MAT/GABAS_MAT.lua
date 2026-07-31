-- AUTHOR: NACHO
-- DATE, TIME WRITEN: 2026/07/30, 09:40, SEMI CLUDY MORNING
-- PLACE OF WRITING: DESAMPARADOS CENTRO, SAN JOSE, SANSCRICIA
-- PURPOSE:
-- GENERAL INFO:
-- SECURITY CONCERNS: NONE
-- LICENSE:
-- DATE, TIME WHEN FIRST STABLE VERSION REACHED:

-- LATER ON I SHALL FINISH THE COBOLESQUE DOCUMENTATION


-- GABAS_MAT.lua
--
-- Every domain lives under modules/, so every require string below is
-- prefixed "modules." -- the same convention GABAS_LINAL.lua itself uses
-- for its own submodules, one level down. Domains that are still just an
-- empty placeholder file (most of them, for now) load fine and return an
-- empty table -- you'll get e.g. GM.Trig == {} until GABAS_TRIG.lua
-- actually has something in it.
--
-- GABAS_LINAL.lua is itself a self-contained package: its own internal
-- requires ("GABAS-LINAL.modules.core", etc.) assume GABAS-LINAL is a
-- direct subfolder of wherever Lua is searching -- exactly true when it's
-- used standalone, but not true once it's nested three levels down inside
-- this master module (GABAS-MAT/modules/GABAS-LINAL/...). Rather than
-- rewriting GABAS_LINAL.lua's internals to know about GABAS-MAT (which
-- would break it for anyone using it standalone), this line adds
-- GABAS-MAT's own modules/ folder as an extra place Lua looks -- so
-- "GABAS-LINAL.modules.core" also resolves under modules/GABAS-LINAL/modules/core.lua.
package.path = package.path .. ";modules/?.lua;modules/?/init.lua"

return {
    Abs_Alg = require("modules.GABAS-ABS-ALG.GABAS_ABS_ALG"),
    Alg = require("modules.GABAS-ALG.GABAS_ALG"),
    Bool_Alg = require("modules.GABAS-BOOL-ALG.GABAS_BOOL_ALG"),
    Calc_One_Var = require("modules.GABAS-CALC-ONE-VAR.GABAS_CALC_ONE_VAR"),
    Calc_Several_Var = require("modules.GABAS-CALC-SEVERAL-VAR.GABAS_CALC_SEVERAL_VAR"),
    Combinatorics = require("modules.GABAS-COMBINATORICS.GABAS_COMBINATORICS"),
    Complex_Anal = require("modules.GABAS-COMPLEX-ANAL.GABAS_COMPLEX_ANAL"),
    Data_Science = require("modules.GABAS-DATA-SCIENCE.GABAS_DATA_SCIENCE"),
    Diff_Eq = require("modules.GABAS-DIFF-EQ.GABAS_DIFF_EQ"),
    Discrete_Math = require("modules.GABAS-DISCRETE-MATH.GABAS_DISCRETE_MATH"),
    Dynamical_Systems = require("modules.GABAS-DYNAMICAL-SYSTEMS.GABAS_DYNAMICAL_SYSTEMS"),
    Encryptation = require("modules.GABAS-ENCRYPTATION.GABAS_ENCRPYTATION"),
    Eng_Math = require("modules.GABAS-ENG-MATH.GABAS_ENG_MATH"),
    Formal_Math = require("modules.GABAS-FORMAL-MATH.GABAS_FORMAL_MATH"),
    Func_Anal = require("modules.GABAS-FUNC-ANAL.GABAS_FUNC_ANAL"),
    Fuzzy_Logic = require("modules.GABAS-FUZZY-LOGIC.GABAS_FUZZY_LOGIC"),
    Game_Thry = require("modules.GABAS-GAME-THRY.GABAS_GAME_THRY"),
    Geom = require("modules.GABAS-GEOM.GABAS_GEOM"),
    Lambda_Calc = require("modules.GABAS-LAMBDA-CALC.GABAS_LAMBDA_CALC"),
    Linal = require("modules.GABAS-LINAL.GABAS_LINAL"),
    Logic = require("modules.GABAS-LOGIC.GABAS_LOGIC"),
    Machine_Learning = require("modules.GABAS-MACHINE-LEARNING.GABAS_MACHINE_LEARNING"),
    Math_Optimization = require("modules.GABAS-MATH-OPTIMIZATION.GABAS_MATH_OPTIMIZATION"),
    Part_Diff_Eqq = require("modules.GABAS-PART-DIFF-EQ.GABAS_PART_DIFF_EQ"),
    Num_Mthd = require("modules.GABAS-NUM-METHD.GABAS_NUM_METHD"),
    Num_Anal = require("modules.GABAS-NUM-ANAL.GABAS_NUM_ANAL"),
    Precalc = require("modules.GABAS-PRECALC.GABAS_PRECALC"),
    Probability = require("modules.GABAS-PROBABILITY.GABAS_PROBABILITY"),
    Prepositional_Calc = require("modules.GABAS-PROPOSITIONAL-CALC.GABAS_PROPOSITONAL_CALC"),
    Real_Anal = require("modules.GABAS-REAL-ANAL.GABAS_REAL_ANAL"),
    Relat_Alg = require("modules.GABAS-RELAT-ALG.GABAS_REL_ALG"),
    Relat_Calculus = require("modules.GABAS-RELAT-CALC.GABAS_RELAT_CALC"),
    Statitstics = require("modules.GABAS-STATISTICS.GABAS_STATISTICS"),
    Stochastic_Calc = require("modules.GABAS-STOCHASTIC-CALC.GABAS_STOCHASTIC_CALC"),
    Tensor_Calc = require("modules.GABAS-TENSOR-CALC.GABAS_TENSOR_CALC"),
    Toplgy = require("modules.GABAS-TOPLGY.GABAS_TOPLGY"),
    Trig = require("modules.GABAS-TRIG.GABAS_TRIG"),
    Vec_Calculus = require("modules.GABAS-VEC-CALC.GABAS_VEC_CALC"),
}


-- B"H.
