/**
 * constants.h
 */
#ifndef __LB_D2Q9_CONSTANTS_H
#define __LB_D2Q9_CONSTANTS_H

#define Lx 64
#define Ly 64

#define tau 0.55
#define o_tau (1.0/tau)
#define o_m_o_tau (1.0 - o_tau)

// ...

#define LB_TYPE FLUIDS

#define EVOLUTION_ALGORITHM TWO_STEP

#endif // __LB_D2Q9_CONSTANTS_H