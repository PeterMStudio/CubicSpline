
#pragma once

#include <QtCore/QList>

#define ErrorGetResult -888.88f


class CubicSplinefromJs {
public:
    enum SplineType{
        Natural = 0, //default
        Quadratic,
        Notaknot,
        Periodic,
    };

private:
    struct ResultData{
        double a;
        double b;
        double c;
        double d;

        double min;
        double max;

        ResultData(double a,double b,double c,double d,double min,double max)
            : a(a),b(b),c(c),d(d),min(min),max(max) { }
    };
public:
    CubicSplinefromJs();
    ~CubicSplinefromJs();

public:
    int setUserPoints(const QList<std::pair<double,double>> &points);

    double getRefY(const double & x);

    void setMinY(const double & minY);
    const double minY() const;

    void setMaxY(const double & maxY);
    const double maxY() const;

private:
    int processPoints( QList<std::pair<double,double>> &points );

    //return val: minX, maxX , minY , maxY
    std::tuple<double,double,double,double>
    getCriticalVal( const QList<std::pair<double,double>> &points);

    //Cubic spline interpolation
    //The function uses the library math.js to ensure high precision results.
    //@param points ,An array of objects with x and y coordinate.
    //@param type The interpolation boundary condition ("quadratic", "notaknot", "periodic", "natural"). "natural" is the default value.
    std::vector<ResultData>  cubicSplineInterpolation( QList<std::pair<double,double>> &points , SplineType type = Natural);

    void outputData( const std::vector<ResultData> & datas);
    // Reduced row echelon form
    // Taken from https://rosettacode.org/wiki/Reduced_row_echelon_form
    // Modified to work with math.js (high float precision).
    std::vector<std::vector<double> >  rref(std::vector<std::vector<double> >  mPoints);

private:
    QList<std::pair<double,double>> m_UserPoints;
    std::vector<ResultData> m_refDatas;
    bool m_bInitSuccess = false;

    bool m_bSetMinY = false;
    bool m_bSetMaxY = false;

    double m_fMinY = 0.0f;
    double m_fMaxY = 255.0f;

};
