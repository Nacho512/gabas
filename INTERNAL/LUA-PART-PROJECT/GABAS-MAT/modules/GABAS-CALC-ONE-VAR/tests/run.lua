local source = debug.getinfo(1, "S").source:sub(2)
local root = source:match("^(.*)/tests/run%.lua$")
local modules_root = root:match("^(.*)/GABAS%-CALC%-ONE%-VAR$")
package.path = modules_root .. "/?.lua;" .. modules_root .. "/?/init.lua;" .. package.path

local tests = {
    "test_calc_one_var",
    "test_quadrature",
    "test_integral",
    "test_inverse_trig",
    "test_integral_infinite",
    "test_integral_principal_value",
    "test_integral_tanh_sinh",
}

for _, name in ipairs(tests) do
    local test = assert(loadfile(root .. "/tests/" .. name .. ".lua"))
    test()
    io.write("PASS ", name, "\n")
end

io.write("PASS ", #tests, " test files\n")
