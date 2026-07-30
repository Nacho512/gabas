-- [FFT]
-- Radix-2 Cooley-Tukey Fast Fourier Transform, O(N log N) instead of the
-- O(N^2) a direct DFT sum would cost. Recursive divide-and-conquer, same
-- shape as Determinant's strassen_mul: split v into its even- and
-- odd-indexed samples, recursively transform each half, then combine with a
-- single "butterfly" pass. Every combine step multiplies by a twiddle
-- factor built as a genuine Complex value (via math.cos/math.sin), so --
-- exactly like the rest of this module -- a real-valued input vector picks
-- up Complex entries "for free" through +/-/* the moment a twiddle factor
-- touches it; FFT never special-cases real vs. Complex input itself.
--
-- Only handles lengths that are an exact power of 2, reusing Core.next_pow2
-- to check that: the recursive halving here has to land on exactly 1 at the
-- bottom, which only happens for a power-of-2 start. This is a real, honest
-- limitation (not implementing Bluestein's algorithm or similar to support
-- arbitrary lengths) -- zero-pad v to the next power of 2 yourself if you
-- need to transform a different length; note that this changes what is
-- actually being transformed (a padded signal's spectrum is not simply "the
-- original spectrum at more points"), so it is left as an explicit choice
-- for the caller to make, not something FFT does silently.
local Core = require("GABAS-LINAL.modules.core")
local Complex = require("GABAS-LINAL.modules.complex")

-- NACHO: THIS MODULE IS INDISPENSABLE FOR THE ALN [ALGEBRA LINEAL] MODULE AND IN ESPECIAL FOR THE TENSOR MODULE.

-- sign = -1 for the forward transform's e^{-i*2*pi*k*n/N} twiddle factors,
-- +1 for the inverse transform's e^{+i*2*pi*k*n/N} ones; IFFT additionally
-- divides every output by N afterwards, which this helper deliberately
-- does not do, since it would be wrong to apply twice across the recursion.
local function fft_recursive(v, sign)
    local n = #v
    if n == 1 then
        return { v[1] }
    end
    local half = n // 2
    local even, odd = {}, {}
    for i = 1, half do
        even[i] = v[2 * i - 1]
        odd[i] = v[2 * i]
    end
    local Feven = fft_recursive(even, sign)
    local Fodd = fft_recursive(odd, sign)
    local result = {}
    for k = 1, half do
        local angle = sign * 2 * math.pi * (k - 1) / n
        local w = Complex.Complex(math.cos(angle), math.sin(angle))
        local term = w * Fodd[k]
        result[k] = Feven[k] + term
        result[k + half] = Feven[k] - term
    end
    return result
end

-- Forward FFT: X_k = sum_n v_n * e^(-i*2*pi*k*n/N). `v` may hold plain
-- numbers, Complex values, or a mix of both.
local function FFT(v)
    assert(type(v) == "table" and #v > 0, "FFT: v must be a non-empty vector.")
    local n = #v
    assert(Core.next_pow2(n) == n,
        "FFT: v's length (" .. n .. ") must be a power of 2 -- the radix-2 Cooley-Tukey " ..
        "algorithm this implements only ever splits its input exactly in half; zero-pad v up " ..
        "to the next power of 2 first if it isn't already one.")
    return fft_recursive(v, -1)
end

-- Inverse FFT: v_n = (1/N) * sum_k V_k * e^(+i*2*pi*k*n/N), recovering the
-- original (generally Complex-valued) sequence a forward FFT(v) produced.
-- Same power-of-2 length restriction as FFT, for the same reason.
local function IFFT(V)
    assert(type(V) == "table" and #V > 0, "IFFT: V must be a non-empty vector.")
    local n = #V
    assert(Core.next_pow2(n) == n,
        "IFFT: V's length (" .. n .. ") must be a power of 2, for the same reason as FFT.")
    local result = fft_recursive(V, 1)
    for i = 1, n do
        result[i] = result[i] / n
    end
    return result
end

return {
    FFT = FFT,
    IFFT = IFFT,
}
