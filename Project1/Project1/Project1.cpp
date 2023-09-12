#include <boost/multiprecision/cpp_int.hpp>
#include <iostream>

typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<4096, 4096, boost::multiprecision::signed_magnitude, boost::multiprecision::unchecked, void> >  int4096_t;

std::pair<int4096_t, int4096_t> EAA_iterative(int4096_t a, int4096_t b) {
    int4096_t x0 = 1, xn = 1;
    int4096_t y0 = 0, yn = 0;
    int4096_t x1 = 0, y1 = 1;
    int4096_t r = a % b, q = a / b;
    while (r > 0) {
        xn = x0 - q * x1;
        yn = y0 - q * y1;
        x0 = x1, y0 = y1;
        x1 = xn, y1 = yn;
        a = b, b = r;
        r = a % b, q = a / b;
    }
    return { x1, y1 };
}

int4096_t modInverse(int4096_t A, int4096_t M)
{
    int4096_t x, y;
    std::tie(x, y) = EAA_iterative(A, M);

    // m is added to handle negative x
    int4096_t res = (x % M + M) % M;
    return res;
}

int4096_t modular_exponentiation(const int4096_t& base, const int4096_t& exponent, const int4096_t& modulus)
{
    int4096_t result = 1;
    int4096_t base_power = base % modulus;

    // Iterate through the bits of the exponent from left to right
    for (int4096_t i = exponent; i > 0; i >>= 1)
    {
        if (i % 2 == 1)
        {
            result = (result * base_power) % modulus;
        }
        base_power = (base_power * base_power) % modulus;
    }

    return result;
}

int main()
{
    using namespace boost::multiprecision;

    int4096_t p ("33478071698956898786044169848212690817704794983713768568912431388982883793878002287614711652531743087737814467999489");
    int4096_t q ("36746043666799590428244633799627952632279158164343087642676032283815739666511279233373417143396810270092798736308917");
    int4096_t ID ("10703014");
    int4096_t n = p * q;
    int4096_t e ("65537");

    // Calculate private key d
    int4096_t phi = (p - 1) * (q - 1);
    int4096_t d = modInverse(e, phi);
    std::cout << "Private key d = " << d << std::endl;
    std::cout << "check: d * e mod phi(n) = " << d * e % phi << std::endl;
    
    // Signing ID directly with d
    int4096_t s = modular_exponentiation(ID, d, n);
    std::cout << "Signiture = " << s << std::endl;

    std::cout << "Validility = " << (modular_exponentiation(s, e, n) == ID) << std::endl;

    // Signing ID with d using CRT
    int4096_t s_alt;

    int4096_t ID_p = ID % p;
    int4096_t ID_q = ID % q;

    int4096_t d_p = d % (p - 1);
    int4096_t d_q = d % (q - 1);

    int4096_t y_p = modular_exponentiation(ID_p, d_p, p);
    int4096_t y_q = modular_exponentiation(ID_q, d_q, q);

    int4096_t c_p = modInverse(q, p);
    int4096_t c_q = modInverse(p, q);

    s_alt = ((q * c_p) % n * y_p + (p * c_q) % n * y_q) % n;

    std::cout << "Signiture (CRT calculation) = " << s_alt << std::endl;

    // check if equal:
    std::cout << "Validility = " << (modular_exponentiation(s_alt, e, n) == ID) << std::endl;

    return 0;
}