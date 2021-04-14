
#include "CubicSplinefromJs.h"
#include <QtCore/QDebug>

CubicSplinefromJs::CubicSplinefromJs() {

}

CubicSplinefromJs::~CubicSplinefromJs() {

}

int CubicSplinefromJs::setUserPoints(const QList<std::pair<double,double>> &points){
    m_bInitSuccess = false;
    m_UserPoints.clear();
    m_UserPoints = points;

    int ret =  processPoints(m_UserPoints);
    if(ret != 0)
        return ret;
    std::vector<ResultData> refDatas = cubicSplineInterpolation(m_UserPoints,Natural);
    if(refDatas.empty())
        return -1;

    m_refDatas = refDatas;
    m_bInitSuccess = true;
    return 0;
}

double CubicSplinefromJs::getRefY(const double & x){
    if(!m_bInitSuccess)
        return ErrorGetResult;

    for(int i=0;i<m_refDatas.size(); ++i){
        if(m_refDatas[i].min <= x && m_refDatas[i].max >= x){
            double y =  m_refDatas[i].a * std::pow(x,3) +
                        m_refDatas[i].b * std::pow(x,2) +
                        m_refDatas[i].c * std::pow(x,1) +
                        m_refDatas[i].d * std::pow(x,0) ;
            if(m_bSetMinY)
                y = y <= m_fMinY ? m_fMinY : y ;
            if(m_bSetMaxY)
                y = y >= m_fMaxY ? m_fMaxY : y ;
            return y;
        }
    }

    return ErrorGetResult;
}

void CubicSplinefromJs::setMinY(const double & minY){
    m_bSetMinY = true;
    m_fMinY = minY;
}

const double CubicSplinefromJs::minY() const{
    return m_fMinY;
}

void CubicSplinefromJs::setMaxY(const double & maxY){
    m_bSetMaxY = true;
    m_fMaxY = maxY;
}

const double CubicSplinefromJs::maxY() const{
    return m_fMaxY;
}

int CubicSplinefromJs::processPoints( QList<std::pair<double,double>> &points ){
    std::sort(points.begin(),points.end(),[=](const std::pair<double,double> & v1, const std::pair<double,double> & v2 ){
        return v1.first < v2.first;
    });


    for(int i=0;i<points.size();++i){
        if( i < points.size()-1  && points[i].first == points[i+1].first ){
            // two points have the same x-value
            // check if the y-value is the same
            if(points[i].second == points[i+1].second ){
                //remove the latter
                points.removeAt(i+1);
            }
            else
            {
                qDebug()<<"Error Points, Same X but Different Y! Error data: x: "<<points[i].first <<" y: "<<points[i].second;
                return -1;
            }
        }
    }

    if(points.size()<2){
        qDebug()<<"Error! Not Enough Points!";
        return -2;
    }

    return 0;

}

std::tuple<double,double,double,double>
CubicSplinefromJs::getCriticalVal( const QList<std::pair<double,double>> &points){
    double minX = points.at(0).first;
    double maxX = points.at(0).first;
    double minY = points.at(0).second;
    double maxY = points.at(0).second;


    for(int i=0 ;i<points.size();++i){
        minX = std::min(minX,points.at(i).first);
        maxX = std::max(maxX,points.at(i).first);
        minY = std::min(minY,points.at(i).second);
        maxY = std::max(maxY,points.at(i).second);
    }

    return std::make_tuple(minX,maxX,minY,maxY);
}

std::vector<CubicSplinefromJs::ResultData> CubicSplinefromJs::cubicSplineInterpolation( QList<std::pair<double,double>> &points , SplineType type){
    int row = 0;
    int solutionIndex = (points.size()-1)*4;

    //Create a dynamic two-dimensional array
    // columns (rows + 1)
    std::vector<std::vector<double> > mPoints( solutionIndex, std::vector<double>(solutionIndex+1) );
    for(int i=0 ; i<solutionIndex; ++i){
        for(int j=0;j<=solutionIndex; ++j){
            mPoints[i][j] = 0.0f;
        }
    }

    // splines through points equations
    for(int i=0 ;i<points.size()-1; ++i,++row){
        std::pair<double,double> point1 = points.at(i);
        std::pair<double,double> point2 = points.at(i+1);
        mPoints[row][i*4+0] = std::pow(point1.first,3);
        mPoints[row][i*4+1] = std::pow(point1.first,2);
        mPoints[row][i*4+2] = point1.first;
        mPoints[row][i*4+3] = 1;
        mPoints[row][solutionIndex] = point1.second;

        mPoints[++row][i*4+0] = std::pow(point2.first,3);
        mPoints[row][i*4+1] = std::pow(point2.first,2);
        mPoints[row][i*4+2] = point2.first;
        mPoints[row][i*4+3] = 1;
        mPoints[row][solutionIndex] = point2.second;
    }

    // first derivative
    for (int functionNr = 0; functionNr < points.size() - 2; functionNr++, row++) {
        std::pair<double,double> point1 = points.at(functionNr+1);
        mPoints[row][functionNr*4+0] = 3*std::pow(point1.first, 2);
        mPoints[row][functionNr*4+1] = 2*point1.first;
        mPoints[row][functionNr*4+2] = 1;
        mPoints[row][functionNr*4+4] = -3*std::pow(point1.first, 2);
        mPoints[row][functionNr*4+5] = -2*point1.first;
        mPoints[row][functionNr*4+6] = -1;
    }

    // second derivative
    for (int functionNr = 0; functionNr < points.size() - 2; functionNr++, row++) {
        std::pair<double,double> point1 = points.at(functionNr+1);
        mPoints[row][functionNr*4+0] = 6*point1.first;
        mPoints[row][functionNr*4+1] = 2;
        mPoints[row][functionNr*4+4] = -6*point1.first;
        mPoints[row][functionNr*4+5] = -2;
    }

    switch (type) {
        case Quadratic:{ // first and last spline quadratic
            mPoints[row++][0] = 1;
            mPoints[row++][solutionIndex-4+0] = 1;
            break;
        }
        case Notaknot :{ // Not-a-knot spline
            mPoints[row][0+0] = 1;
            mPoints[row++][0+4] = -1;
            mPoints[row][solutionIndex-8+0] = 1;
            mPoints[row][solutionIndex-4+0] = -1;
            break;
        }
        case Periodic:{
            // first derivative of first and last point equal
            mPoints[row][0] = 3* std::pow( points[0].first,2);
            mPoints[row][1] = 2* points[0].first;
            mPoints[row][2] = 1;
            mPoints[row][solutionIndex-4+0] =  -3* std::pow(points[points.size()-1].first,2) ;
            mPoints[row][solutionIndex-4+1] =  -2* points[points.size()-1].first;
            mPoints[row++][solutionIndex-4+2] = -1;

            // second derivative of first and last point equal
            mPoints[row][0] = 6*points[0].first;
            mPoints[row][1] = 2;
            mPoints[row][solutionIndex-4+0] = -6 * points[points.size()-1].first;
            mPoints[row][solutionIndex-4+1] = -2;
            break;
        }
        case Natural:{
            mPoints[row][0+0] = 6*points[0].first;
            mPoints[row++][0+1] = 2;
            mPoints[row][solutionIndex-4+0] = 6*points[points.size()-1].first;
            mPoints[row][solutionIndex-4+1] = 2;
            break;
        }
    }

    auto reducedRowEchelonForm =  rref(mPoints);

    std::vector<ResultData> vectorData;
    std::vector<double> coefficients;
    for(int i=0;i<reducedRowEchelonForm.size();++i){
        coefficients.push_back( reducedRowEchelonForm[i][reducedRowEchelonForm[i].size()-1]);
    }

    if(reducedRowEchelonForm.empty())
        return vectorData;


    for(int i=0;i<coefficients.size();i+=4){
        ResultData data(coefficients[i],  coefficients[i+1],
                        coefficients[i+2],coefficients[i+3],
                        points[i/4].first,points[i/4+1].first );
        vectorData.push_back(data);
    }

    return vectorData;
}

void CubicSplinefromJs::outputData(const std::vector<ResultData> & datas){
    for(int i=0;i<datas.size();++i){
        qDebug()<<QString("range from %1 to %2").arg(datas[i].min).arg(datas[i].max);
        qDebug()<<QString( "%1 * x^3 + (%2) * x^2  +(%3) * x^1 + (%4)*x^0 " )
                    .arg(datas[i].a)
                    .arg(datas[i].b)
                    .arg(datas[i].c)
                    .arg(datas[i].d)
                    ;
    }
}

std::vector<std::vector<double> >  CubicSplinefromJs::rref(std::vector<std::vector<double> >  Points){

    std::vector<std::vector<double> > mPoints = Points;
    int lead = 0;
    for(int r = 0; r<mPoints.size() ; ++r){
        if(mPoints[0].size()<= lead)
            return std::vector<std::vector<double> >{};
        int i = r;
        while (mPoints[i][lead] == 0){
            ++i;
            if(mPoints.size() == i){
                i = r;
                ++lead;
                if(mPoints[0].size() == lead)
                    return std::vector<std::vector<double> >{} ;
            }
        }

        auto tmp = mPoints[i];
        mPoints[i] = mPoints[r];
        mPoints[r] = tmp;

        double val = mPoints[r][lead];
        for(int j=0;j<mPoints[0].size();++j){
            mPoints[r][j] = mPoints[r][j] / val;
        }

        for(int i=0; i<mPoints.size() ; ++i){
            if(i == r)
                continue;
            val = mPoints[i][lead];
            for(int j=0;j<mPoints[0].size() ; ++j){
                mPoints[i][j] = mPoints[i][j] - val * mPoints[r][j];
            }
        }
        lead++;
    }
    return mPoints;
}
