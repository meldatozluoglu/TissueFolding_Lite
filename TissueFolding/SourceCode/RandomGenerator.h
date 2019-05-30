/*
 * RandomGenerator.h
 *
 *  Created on: 5 Jan 2016
 *      Author: melda
 */

#ifndef RANDOMGENERATOR_H_
#define RANDOMGENERATOR_H_


#include <boost/random.hpp>
#include <iostream>
#include <ctime>
#include <vector>

/**
 * The random number generator.
 */
class RandomGenerator
{

public:

    /**
     * Default constructor.
     */
    RandomGenerator()
    : m_rnGen (static_cast< unsigned int > (std::time(0)))
    {
        this->Seed();
    }

    /**
     * Get a given number of normally-distributed random variables.
     */
    std::vector <double> getNormRV( double mean, double var, unsigned int n )
    {
        std::vector<double> res;

        // create the normal distribution for rotational
        // velocities
        boost::normal_distribution< double > distrib ( mean, var );

        // create the var generator.
        boost::variate_generator< boost::mt19937&,
            boost::normal_distribution<double> > rvGen (
            m_rnGen, distrib );

        for ( unsigned int i = 0; i < n; ++i )
        {
            // prevent occasional huge values
            // define a boundary that is 1000 times deviation
            double _rv = rvGen();
            double _bound = sqrt(var) * 100;

            // if out of boundary, pick another
            while ( _rv < -_bound || _rv > _bound )
            {
                _rv = rvGen();
            }

            res.push_back( _rv );
        }

        return res;
    }

    /**
     * Get a normally-distributed random variable.
     */
    double getNormRV( double mean, double var )
    {
        // create the normal distribution for rotational
        // velocities
        boost::normal_distribution< double > distrib ( mean, var );

        // create the var generator.
        boost::variate_generator< boost::mt19937&,
            boost::normal_distribution<double> > rvGen (
            m_rnGen, distrib );

        // prevent occasional huge values
        // define a boundary that is 100 times deviation
        double _rv = rvGen();
        double _bound = sqrt(var) * 100;

        // if out of boundary, pick another
        while ( _rv < -_bound || _rv > _bound )
        {
            _rv = rvGen();
        }

        return _rv;
    }

    /**
     * Get a uniformly distributed RV in [m,n] (m,n are int).
     */
    int getUniformRV( int min, int max )
    {
        // create the distribution
        boost::uniform_int<> distro ( min, max );

        // create the rv gen
        boost::variate_generator< boost::mt19937&,
            boost::uniform_int<> > rvGen ( m_rnGen, distro );

        return rvGen();
    }

    /**
     * Get a uniformly distributed RV in [m,n] (m,n are int).
     */
    double getUniformRV( double min, double max )
    {
        // create the distribution
        boost::uniform_real<> distro ( min, max );

        // create the rv gen
        boost::variate_generator< boost::mt19937&,
            boost::uniform_real<> > rvGen ( m_rnGen, distro );

        return rvGen();
    }

    /**
     * Get a randomised string of a certain length.
     */
    std::string getRandomStr( unsigned int length )
    {
        std::string _res;

        for ( unsigned int i = 0; i < length; ++i )
        {
            int _c = getUniformRV( 97, 122 );
            _res.push_back( (char)_c );
        }

        return _res;
    }

    /**
     * Seed the generator.
     */
    void Seed()
    {
        m_rnGen.seed();
    }

    /**
     * The random generator access function.
     */
    static RandomGenerator& Obj()
    {
        if ( !m_pObj )
        {
            m_pObj = new RandomGenerator();
            m_pObj->m_rnGen.seed( std::time(0) );
        }
        return *m_pObj;
    }

private:

    /**
     * The boost Mersenne-Twister random number generator.
     */
    boost::mt19937 m_rnGen;

    /**
     * The global singleton.
     */
    static RandomGenerator* m_pObj;
};



#endif /* RANDOMGENERATOR_H_ */
