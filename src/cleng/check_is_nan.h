#include "cleng.h"

using Fp_info = numeric_limits<double>;
inline auto is_ieee754_nan( double const x ) -> bool {
    static constexpr bool   is_claimed_ieee754  = Fp_info::is_iec559;
    static constexpr int    n_bits_per_byte     = CHAR_BIT;
    using Byte = unsigned char;

    static_assert( is_claimed_ieee754, "!" );
    static_assert( n_bits_per_byte == 8, "!" );
    static_assert( sizeof( x ) == sizeof( uint64_t ), "!" );

#ifdef _MSC_VER
    uint64_t const bits = reinterpret_cast<uint64_t const&>( x );
#else
    Byte bytes[sizeof(x)];
    memcpy( bytes, &x, sizeof( x ) );
    uint64_t int_value;
    memcpy( &int_value, bytes, sizeof( x ) );
    uint64_t const& bits = int_value;
#endif

    static constexpr uint64_t   sign_mask       = 0x8000000000000000;
    static constexpr uint64_t   exp_mask        = 0x7FF0000000000000;
    static constexpr uint64_t   mantissa_mask   = 0x000FFFFFFFFFFFFF;

    (void) sign_mask;
    return (bits & exp_mask) == exp_mask and (bits & mantissa_mask) != 0;
}
