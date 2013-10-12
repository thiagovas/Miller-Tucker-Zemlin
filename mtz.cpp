#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <utility>
#include <cmath>
#include <set>
#include <map>
#include "ilcplex/ilocplex.h"
using namespace std;
ILOSTLBEGIN

void PrintVector(vector<int> vec);

/*
	Implementação do Miller-Tucker-Zemlin
	Feito por Thiago Vieira
*/

// t representa o número de retornos para a cidade origem zero.
#define T 1

int main()
{
	ios::sync_with_stdio(false);
	IloEnv env;
	IloModel model(env);
	
	int n, peso;
	
	bool AddTRestriction = true;
	
	try{
		cin >> n;
		IloArray<IloIntVarArray> x(env, n);
		IloArray<IloIntArray> d(env, n);
		
		for(int i = 0; i < n; i++)
		{
			d[i] = IloIntArray(env);
			for(int j = 0; j < n; j++)
			{
				cin >> peso;
				d[i].add(peso);
			}
		}
		
		for(int i = 0; i < n; i++)
			x[i] = IloIntVarArray(env, n, 0, 1);
		
		
		if(AddTRestriction)
		{
			/* BEGIN Restrição do número de retornos para a cidade origem. */
			IloExpr exprTempPrim(env); // Restrição sobre o número de arestas incidentes sobre o vértice origem.
			
			for(int i = 1; i < n; i++)
				exprTempPrim += x[i][0];
			
			model.add(exprTempPrim == T);
			exprTempPrim.end();
			/* END Restrição do número de retornos para a cidade origem. */
		}
		
		for(int i = 1; i < n; i++)
		{
			IloExpr exprTempPrim(env); // Restrição sobre o número de arestas incidentes.
			IloExpr exprTempSecond(env); // Restrição sobre o número de arestas saindo do vértice.
			
			for(int j = 0; j < n; j++)
			{
				if(i == j) continue;
				exprTempPrim += x[i][j];
				exprTempSecond += x[j][i];
			}
			
			model.add(exprTempPrim == 1);
			model.add(exprTempSecond == 1);
			exprTempPrim.end();
			exprTempSecond.end();
		}
		
		
		/* BEGIN - LEVEL VARIABLES RESTRICTIONS */
		IloIntVar v0 = 0;
		IloIntVarArray levelu(env, n, 0, n);
		levelu[0] = v0;
		
		for(int i = 1; i < n; i++)
		{
			for(int j = 1; j < n; j++)
			{
				if(i == j) continue;
				
				model.add(levelu[i]-levelu[j] + n*x[i][j] <= n-1);
			}
		}
		/* END - LEVEL VARIABLES */
		
		
		/* BEGIN - OBJECTIVE FUNCTION */
		IloExpr obj(env);
		for(int i = 0; i < n; i++)
		{
			obj += IloScalProd(d[i], x[i]);
			obj -= d[i][i]*x[i][i];
		}
		model.add(IloMinimize(env, obj));
		obj.end();
		/* END - OBJECTIVE FUNCTION */
		
		cout << endl;
		IloCplex cplex(model);
		
		if(cplex.solve())
		{
			env.out() << "\nSolution Status = " << cplex.getStatus() << endl;
			env.out() << "Solution value = " << cplex.getObjValue() << endl;
			
			// Printando os valores das variáveis.
			for(int i = 0; i < n; i++)
			{
				for(int j = 0; j < n; j++)
					if(i == j)
					{
						cout << "0 ";
						continue;
					}
					else cout << cplex.getValue(x[i][j]) << " ";
				cout << endl;
			}
		}
		else
		{
			env.out() << "Solution Status = " << cplex.getStatus() << endl;
			cout << "Deu pau, nao resolveu!\n";
		}
	}
	catch(IloException &e)
	{
		cerr << "OH: " << e << endl;
	}
	catch(...)
	{
		cerr << "Alguma outra coisa aconteceu!\n";
	}
	
	env.end();
	return 0;
}

/* Método para imprimir qualquer vector de ints. */
void PrintVector(vector<int> vec)
{
	if(vec.empty()) return;
	
	cout << *vec.begin();
	for(vector<int>::iterator it = ++vec.begin(); it != vec.end(); it++)
		cout << " " << *it;
	
	cout << endl;
}
