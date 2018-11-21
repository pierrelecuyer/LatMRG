#include "latmrg/mrgtypes/MixmaxMMRG.h"
#include "latticetester/Util.h"

using namespace std;
using namespace LatticeTester;

namespace LatMRG
{

//===========================================================================

MixmaxMMRG::MixmaxMMRG(MScal modulus, int & N, MScal & s, MScal & m, MScal & b)
{
	m_modulus = modulus;
	m_order = N;
	m_parameter1 = s;
	m_parameter2 = m;
	m_parameter3 = b;
	buildMatrix(m_modulus, m_order, m_parameter1, m_parameter2, m_parameter3);

}

//===========================================================================

MixmaxMMRG::MixmaxMMRG(MScal modulus, int & N, MScal & s, MScal & m)
{
	m_modulus = modulus;
	m_order = N;
	m_parameter1 = s;
	m_parameter2 = m;
	m_parameter3 = 2 - 2*m;
	buildMatrix(m_modulus, m_order, m_parameter1, m_parameter2, m_parameter3);

}

//===========================================================================

MixmaxMMRG::MixmaxMMRG(MScal modulus, int & N, MScal & s)
{
	m_modulus = modulus;
	m_order = N;
	m_parameter1 = s;
	m_parameter2 = 1;
	m_parameter3 = 0; // =2-2*m
	buildMatrix(m_modulus, m_order, m_parameter1, m_parameter2, m_parameter3);

}

//===========================================================================

MixmaxMMRG::~MixmaxMMRG()
{}

//===========================================================================

void MixmaxMMRG::buildMatrix(MScal modulus, int & N, MScal & s, MScal & m, MScal & b)
{
	m_A.kill();
    m_A.resize(N,N);
    for (int j = 1; j < N; j ++) {
        for (int i = j+1; i < N; i++) {
            m_A[i][j] = (i-j+2) * m + b;
        	Modulo(m_A[i][j], modulus, m_A[i][j]);
        }
    }
    for (int i = 0; i < N; i ++)
        m_A[i][0] = 1;
    for (int i = 1; i < N; i++)
        m_A[i][i] = 2;
    for (int i = 0; i < N; i ++) {
        for (int j = i+1; j < N; j ++)
            m_A[i][j] = 1;
    }
    m_A[2][1] += s;
    Modulo(m_A[2][1], modulus, m_A[2][1]);
}

//===========================================================================
   
} // end namespace
