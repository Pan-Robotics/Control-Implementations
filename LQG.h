#ifndef LQG_H
#define LQG_H

#include <iostream>
#include <Eigen/Dense>

class LQG
{
private:
    Eigen::MatrixXd Adt;
    Eigen::MatrixXd Bdt;
    Eigen::MatrixXd Cdt;
    Eigen::MatrixXd Ddt;
    Eigen::MatrixXd Kdt;
    Eigen::MatrixXd Kidt;
    Eigen::MatrixXd Ldt;
    Eigen::MatrixXd U_e;

    Eigen::MatrixXd U;
    Eigen::MatrixXd prevU;

    Eigen::MatrixXd Xest;
    Eigen::MatrixXd prevXest;

    Eigen::MatrixXd Xe;
    Eigen::MatrixXd prevXe;

    Eigen::MatrixXd Ref;
    Eigen::MatrixXd Y;
    Eigen::MatrixXd e;

public:
    // Constructor
    LQG(const Eigen::MatrixXd &Adt, const Eigen::MatrixXd &Bdt, const Eigen::MatrixXd &Cdt, const Eigen::MatrixXd &Ddt,
        const Eigen::MatrixXd &Kdt, const Eigen::MatrixXd &Kidt, const Eigen::MatrixXd &Ldt, const Eigen::MatrixXd &U_e)
    {
        this->Adt = Adt;
        this->Bdt = Bdt;
        this->Cdt = Cdt;
        this->Ddt = Ddt;
        this->Kdt = Kdt;
        this->Kidt = Kidt;
        this->Ldt = Ldt;
        this->U_e = U_e;

        U = Eigen::MatrixXd::Zero(Bdt.rows(), 1);
        prevU = Eigen::MatrixXd::Zero(Bdt.rows(), 1);

        Xest = Eigen::MatrixXd::Zero(Adt.rows(), 1);
        prevXest = Eigen::MatrixXd::Zero(Adt.rows(), 1);

        Xe = Eigen::MatrixXd::Zero(Cdt.rows(), 1);
        prevXe = Eigen::MatrixXd::Zero(Cdt.rows(), 1);

        Ref = Eigen::MatrixXd::Zero(Cdt.rows(), 1);
        Y = Eigen::MatrixXd::Zero(Cdt.rows(), 1);
        e = Eigen::MatrixXd::Zero(Cdt.rows(), 1);
    }

    // Calculate method
    Eigen::MatrixXd calculate(const Eigen::MatrixXd &prevU, const Eigen::MatrixXd &Y, const Eigen::MatrixXd &Ref, bool Linear)
    {
        this->Y = Y;
        this->prevU = prevU;
        this->Ref = Ref;

        if (Linear)
        {
            Xest = Adt * prevXest + Bdt * prevU; // KF Linear Prediction
            e = Y - Xest(Eigen::seq(0, 10, 2), Eigen::all);
            Xest += Ldt * e;

            Xe = prevXe + (Ref - Xest(Eigen::seq(0, 10, 2), Eigen::all)); // Integrator
            U = -(Kdt * Xest) - (Kidt * Xe);
        }
        else
        {
            Xest = Adt * prevXest + Bdt * (prevU - U_e); // KF Non-Linear Prediction
            e = Y - Xest(Eigen::seq(0, 10, 2), Eigen::all);
            Xest += Ldt * e;

            Xe = prevXe + (Ref - Xest(Eigen::seq(0, 10, 2), Eigen::all)); // Integrator
            U = (U_e - (Kdt * Xest) - (Kidt * Xe)).cwiseMax(0).cwiseMin(800);
        }

        prevXe = Xe;
        prevXest = Xest;

        return U;
    }
};

#endif // LQG_H
