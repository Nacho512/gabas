local Core = require("GABAS-LINAL.modules.core")

local function Assert_Close(actual, expected, tolerance, label)
    tolerance = tolerance or Core.DEFAULT_ABS_TOL
    Core.assert_nonneg_number(tolerance, "Assert_Close: tolerance")
    Core.assert_finite_scalar(actual, "Assert_Close: actual")
    Core.assert_finite_scalar(expected, "Assert_Close: expected")
    local difference = Core.cabs(actual - expected)
    assert(difference <= tolerance,
        (label or "values") .. " differ by " .. difference .. " (tolerance " .. tolerance .. ").")
end

local function Assert_Vector_Close(actual, expected, tolerance, label)
    local n = Core.assert_vector(actual, "Assert_Vector_Close: actual")
    Core.assert_vector(expected, "Assert_Vector_Close: expected", {length = n})
    for i = 1, n do
        Assert_Close(actual[i], expected[i], tolerance, (label or "vectors") .. " at [" .. i .. "]")
    end
end

local function Assert_Matrix_Close(actual, expected, tolerance, label)
    local rows, cols = Core.assert_matrix(actual, "Assert_Matrix_Close: actual")
    local expected_rows, expected_cols = Core.assert_matrix(expected, "Assert_Matrix_Close: expected")
    assert(rows == expected_rows and cols == expected_cols,
        "Assert_Matrix_Close: matrices must have matching dimensions.")
    for i = 1, rows do
        for j = 1, cols do
            Assert_Close(actual[i][j], expected[i][j], tolerance,
                (label or "matrices") .. " at [" .. i .. "][" .. j .. "]")
        end
    end
end

local function Assert_Error(fn, pattern)
    assert(type(fn) == "function", "Assert_Error: fn must be a function.")
    local ok, message = pcall(fn)
    assert(not ok, "Assert_Error: expected an error, but the call succeeded.")
    if pattern ~= nil then
        assert(type(pattern) == "string", "Assert_Error: pattern must be a string.")
        assert(tostring(message):find(pattern, 1, true),
            "Assert_Error: error did not contain expected text '" .. pattern .. "': " .. tostring(message))
    end
end

return {
    Assert_close = Assert_Close,
    Assert_vector_close = Assert_Vector_Close,
    Assert_matrix_close = Assert_Matrix_Close,
    Assert_error = Assert_Error,
}
