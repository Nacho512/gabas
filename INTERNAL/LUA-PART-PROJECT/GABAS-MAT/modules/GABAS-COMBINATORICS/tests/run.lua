local source = debug.getinfo(1, "S").source:sub(2)
local root = source:match("^(.*)/tests/run%.lua$")
local modules_root = root:match("^(.*)/GABAS%-COMBINATORICS$")
package.path = modules_root .. "/?.lua;" .. modules_root .. "/?/init.lua;" .. package.path

local tests = {
    "test_combinatorics",
}

for _, name in ipairs(tests) do
    local test = assert(loadfile(root .. "/tests/" .. name .. ".lua"))
    test()
    io.write("PASS ", name, "\n")
end

io.write("PASS ", #tests, " test files\n")
