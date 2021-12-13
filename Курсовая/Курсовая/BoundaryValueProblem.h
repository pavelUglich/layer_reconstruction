#include<vector>
#include<complex>


using namespace std;

#pragma once


/**
 * \brief ������� ������ ��� ����������� �������� ���� 
 * (������� ���������� ������������)
 */
class BoundaryValueProblem {
protected:
	
	/**
	 * \brief ���������� ������ ��������
	 */
	double innerRadius;

	/**
	 * \brief ������� ��������������� ������ ����
	 */
	vector<complex<double>> CauchyProblemSolution;

	/**
	 * \brief ����������� ��������
	 */
	double eps;

	/**
	 * \brief �������� �������������� �����
	 */
	//complex<double> alpha;

	/**
	 * \brief �������� �����
	 */
	double kappa;


	
public:

	//vector<complex<double>> Shoot(double a, double b, vector<complex<double>> & InitialValue) const = 0;

	virtual vector<complex<double>> CauchyProblem(double, double, vector<complex<double>>&) const = 0;

	/**
	 * \brief ������� ��������������� ������ ����
	 * \param y 
	 * \param size 
	 * \return 
	 */
	virtual vector<complex<double>> CauchyProblemSolutions(double y,
		size_t size = 4) const = 0;

	vector<vector<complex<double>>> InnerCauchyProblemSolution;	
	BoundaryValueProblem();

	/**
	 * \brief ����������� 
	 * \param innerRadius ���������� ������ ��������
	 * \param kappa �������� �����
	 * \param alpha �������� �������������� �����
	 * \param eps ����������� ����������
	 */
	BoundaryValueProblem(double innerRadius, double kappa, double eps = 0.1e-5);
	/**
	 * \brief ����������
	 */
	virtual ~BoundaryValueProblem();


	/**
	 * \brief ����� �������
	 * \return ������
	 */
	complex<double> NewtonMethod();


	/**
	 * \brief 
	 * \return 
	 */
	vector<complex<double>> GetCauchyProblemSolution() const { return CauchyProblemSolution; };
	
	/**
	 * \brief ����������� ���������������� ���������
	 * \param size ������ ������� ���������������� ���������
	 * \return �����������
	 */
	complex<double> IntegrandDenumerator(size_t size = 2) const;


	/**
	 * \brief �������� ���� �� ������� ����������� ��������
	 * \return �������� ���� �� ������� ����������� ��������
	 */
	virtual complex<double> waveField() const;

	/**
	 * \brief �������� �������� ��������� �����.
	 * \param _alpha 
	 */
	void SetAlpha(complex<double> _alpha);
	
	/**
	 * \brief �������, ��������� �� ���������� ����������
	 * \param points ���������� �����
	 * \param size ���������� ��������� �������
	 */
	virtual void CreateTheInnerSolution(size_t points, size_t size = 4) = 0;

	/**
	 * \brief ������
	 * \return �������� �����
	 */
	double getKappa() const;
	double getq() const;
};

