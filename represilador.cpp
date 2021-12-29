/*
 * @file repressilator.cpp
 * @brief Silmulación del Represilador (Elowitz & Leibler, 2000) en C++
 *
 * Copyright 2021 Juan C. Arboleda R. <juan.arboleda2@udea.edu.co>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <iostream>
#include <math.h>
#include <boost/numeric/odeint.hpp>

namespace odeint = boost::numeric::odeint;

typedef std::vector<double> state_type;

void repressilator(const state_type& x, state_type& dxdt, const double t){
    const double a = 216.4;
    const double a0 = 0.2164;
    const double n = 2.0; // Coeficiente de Hill
    const double b = 0.2; // Deicaimiento proteína/decaimiento mRNA

    const double& m_lacI = x[0];
    const double& p_lacI = x[1];
    const double& m_tetR = x[2];
    const double& p_tetR = x[3];
    const double& m_cI = x[4];
    const double& p_cI = x[5];

    // Ecuación diferencial para el mRNA LacI
    double& dmlacI_dt = dxdt[0] = -m_lacI + a/(1.0 + pow(p_cI, n) ) + a0 ;
    // Ecuación diferencial para la proteína LacI
    double& pdlacI_dt = dxdt[1] = -b*( p_lacI - m_lacI ) ;

    // Ecuación diferencial para el mRNA TetR
    double& dmtetR_dt = dxdt[2] = -m_tetR + a/(1.0 + pow(p_lacI, n) ) + a0 ;
    // Ecuación diferencial para la proteína TetR
    double& dptetR_dt = dxdt[3] = -b*( p_tetR - m_tetR ) ;

    // Ecuación diferencial para el mRNA cI
    double& dmcI_dt = dxdt[4] = -m_cI + a/(1.0 + pow(p_tetR, n) ) + a0 ;
    // Ecuación diferencial para la proteína TetR
    double& dpcI_dt = dxdt[5] = -b*( p_cI - m_cI ) ;
}

// La siguiente función imprime en pantalla el estado del sistema en el tiempo t
void write_out(const state_type& x, const double& t){
    std::cout << t << '\t' << x[1] << '\t' << x[3] << '\t' << x[5] << '\n';
}

// Esta función simplemente imprimirá en pantalla el encabezado de la tabla
// de datos
void write_header(){
    std::cout << "Tiempo\t" << "LacI\t" << "TetR\t" << "cI\n";
}

int main(){
    /* En primer lugar creamos el vector con el estado inicial del sistema.
    * Iniciaremos con 20 unidades del mRNA correspondiente a TetR.
    */
    state_type x0{0, 0, 20, 0, 0, 0};

    // Escribimos el encabezado de la tabla
    write_header();

    // Integramos las ecuaciones diferenciales usando la libería 'odeint'
    odeint::integrate(repressilator,
                      x0,
                      0.0, 500.0, 0.1, // Integraremos desde t=0 hasta t=500, el tamaño de paso inicial será 0.1
                      write_out);
}
