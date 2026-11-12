#include"user_funs.h"

matrix ff0T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego
{
	matrix y;												// y zawiera warto�� funkcji celu
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);		// ud1 zawiera wsp�rz�dne szukanego optimum
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego
{
	matrix y;												// y zawiera warto�� funkcji celu
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
		MT = matrix(2, new double[2] { m2d(x), 0.5 });		// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	matrix* Y = solve_ode(df0, 0, 0.1, 20, Y0, ud1, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	int n = get_len(Y[0]);									// d�ugo�� rozwi�zania
	double teta_max = Y[1](0, 0);							// szukamy maksymalnego wychylenia wahad�a
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));							// warto�� funkcji celu (ud1 to za�o�one maksymalne wychylenie)
	Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double m = 1, l = 0.5, b = 0.5, g = 9.81;				// definiujemy parametry modelu
	double I = m * pow(l, 2);
	dY(0) = Y(1);																// pochodna z po�o�enia to pr�dko��
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;	// pochodna z pr�dko�ci to przyspieszenie
	return dY;
}





matrix ff1R(matrix x, matrix ud1, matrix ud2){
	matrix y;
	matrix Y0(3, new double[3]{5,1,20});
	matrix *Y = solve_ode(df1,0,1,2000,Y0,ud1,x);
	int n = get_len(Y[0]);
	double max_temp = Y[1](0,2);
	for(int i=1; i<n; i++){
		if(max_temp < Y[1](i,2))
			max_temp = Y[1](i,2);
	}
	y = abs(max_temp - 50);
	delete[]Y;
	return y;
}


matrix ff1T(matrix x, matrix ud1 , matrix ud2){
	matrix y;
	y = -cos(0.1 * x())*exp(-pow(0.1*x()-2*M_PI, 2))+0.002*pow(0.1*x(),2);
	return y;
}


matrix ff2T(matrix x, matrix ud1, matrix ud2){
	//x wektor zmiennych decyzyjnych:
	//x(0) = x1
	//x(1)=x2
	//
	//ud1 ud2 dodatkowe argumenty (nie wykorzystujemy ich w tej funkcji)
	
matrix y; //zmienna przechowujaca wartosc funkcji celu 

// funkcja testowa z konspektu:
//f(x1,x2) = x1^2 + x2^2 -cos(2.5 * pi *x1) - cos(2.5 * pi *x2) +2
//M_PI stala ozn licbe pi
//funkcja ma wiele min lokalnych
//mis globalne:
//x1 = 0
//x2 = 0
//f(0,0)=0
y = pow(x(0), 2)+ pow(x(1), 2) - cos(2.5 * M_PI * x(0)) - cos(2.5 *M_PI *x (1)) + 2;
return y;
}



matrix ff2R(matrix x, matrix ud1, matrix ud2){
	//x - wektor wsp regulatora 
	//x(0) = k1 wspolczynnik zwiazany z bledem polozenia 
	//x(1) = k2 wspolczynnik zwiazany z bledem predkosci
	//zadaniem funkcji jest oblixznie wartosci funkcji jakosci sterpwania Q(k1,k2)

	matrix y; //zmienna przechowujaca wartosc funkcji celu

	//Y0 wektor warunkow poczatkowych o wymiarach 2x1
	//Y0(0) = 0 poczatkowe polozenie katowe w radianach
	//Y0(1) = 0 pocztakowa predkosc katowa rad/s
	//konstruktor matrix(1,2) wypelnia macierz zerami 
	//Y_ref wektor wartosci zadanych:
	//Y_ref(0) = 3.14 zadany kat obrotu ramienia w radianach
	//Y_ref(1)=0 - zadana predkosc katowa w rad/s

	matrix Y0(2,1), Y_ref(2, new double[2]{3.14, 0});

	//rozwiazujemy rownania rozniczkowe ruchu ramienia gdzie 
	/*
	df2 - funkcja opisujaca model ramienia robota
	0 czas poczatkowy w s
	0.1 krok czasowy w s
	100  czas koncowy w s
	Y0 - warunki poczatkowe 
	Y_ref zadane polozenie i predkosc 
	x - wsp regulatora k1 i k2 
	Y[0] - wektor kolejnych chwil czasu \Y[1] macierz polozenia i predkosi ramienia

	kazde wywolanie ff2r powoduje wykonanie calej symulacji ruchu ramienia robota
	*/

	matrix* Y =solve_ode(df2, 0, 0.1, 100, Y0, Y_ref, x);
	/*
	pobieramy liczbe punktow symulacji \\dla przedzialu 0 - 100 s i kroku 0.1s otrzymujemy 1001 punktow 
	*/
	int n = get_len(Y[0]);
	y = 0; // poczatkowa wartosc funkcji jakosci 
	//przechodzimy przez wszystkie punkty symulacji 
	for (int i = 0; i <n; i++){
		/*obliczamy funkcje jakosci
		Q - calka z:
		10*(alpha_ref - alpha)^2 + (omega_ref - omega)^2 + M(t)^2

		peirwszy skladnik: kara za blad polozenia katowego 
		drugi skladnik: kara za blad predkosci katowej
		trzeci skladnik: kara za duzy moment sterujacy

		Y[1](i,0) - polozenie ramienia w chwili i
		Y[1](i,1) polozenie ramienia w chwili i 

		moment sterujacy:

		M = k1*(alpha_ref - alpha) + k2*(omega_Ref - omega)
		*/

		y = y+10 *pow(Y_ref(0) - Y[1](i,0),2) + pow(Y_ref(1) - Y[1](i,1),2) + pow(x(0) * (Y_ref(0) - Y[1](i,0)) + x(1) * (Y_ref(1) - Y[1](i,1)),2);
	}

	/*
	mnozymy sume przez krok czasowy dt = 0.1
	w ten sposob przyblizamy calke za pomoca metody prostokatow */

	y = y * 0.1;
	delete[] Y;
	/*
	zwalniamy pamiec zajmowana przez wyniki symulacji, zwracamy obliczona wartosc funkcji jakosci 
	
	algorytmy optymalizacji beda probowaly znalezc takie k1 i k2 aby y bylo jak najmnijejsze
	*/

	return y;
}


matrix df2(double t, matrix Y, matrix ud1, matrix ud2){
	/*
	t - aktualny czas symulacji 
	Y(0)- polozenie katowe w radianach (alpha)
	Y(1)- predkosc katowa w rad/s (omega)

	ud1 wartosci zadane:
	ud1(0) alpha_ref - zadany kat 
	ud1(1) omega_ref - zadana predkosc

	ud2 - wspolczynniki regulatora 
	ud2(0) = k1
	ud2(1) = k2 

	parametry fizyczne modelu:
	mr = 1kg - masa ramienia
	mc = 5kg - masa ciezarka 
	l = 2m - dlugosc ramienia
	b = 0.25 - wsp tarcia
	*/
double mr = 1, mc = 5, l = 2, b = 0.25;

/*
moment bezwladnosci ramienia z ciezarkiem:
I = (mr*l^2)/3 + mc*l^2

dla przyjetych parametrow:
I = (1*4)/3 + 5*4
I = 21.3333... kg*m^2
*/
double I = mr * l * l / 3 + mc * l * l;
//tworzymy wektor pochodnych o wymiarach 2x1
//dY(0) -  pochodna polozenia katowego
//dY(1) - pochodna predkosci katowej
matrix dY(2,1);
/*pierwsze rownanie:
d(alpha)/dt + b*omega = M(t)
po przeksztalceniu :
d(omega)/dt = (M(t) - b *omega)

moment sterujacy:  M(t) = k1*(alpha_ref - alpha) + k2*(omega_ref - omega)

cale wyrazeniei oblicza przyspieszenie katowe:
*/
dY(0) = Y(1);
dY(1) = (ud2(0) * (ud1(0) - Y(0)) + ud2(1) * (ud1(1) - Y(1)) - b * Y(1))/I;
//zwracamy wrektor pochodnych ktory zostanie wykorzystany przez solve_ode() do numeryxznego rozwiazania rownan ruchu 

return dY;
}


matrix df1(double t, matrix Y, matrix ud1, matrix ud2){
double a = 0.98, b = 0.63, g = 9.81;
double PA = 1.5, TA = 90, PB = 1;
double DB = 0.00365665, Fin=0.01, Tin=20;
matrix dY (3,1);
double FAout = Y(0)>0?a*b*m2d(ud2)*sqrt((2*g*Y(0))/PA):0;
double FBout = Y(1)>0?a*b*DB*sqrt(2*g*Y(1)/PB):0;
dY(0)= -FAout;
dY(1)= FAout + Fin -FBout;
dY(2)= Fin/Y(1)*(Tin -Y(2))+FAout/Y(1)*(TA - Y(2));
return dY;
}

matrix ff3Ta(matrix x, matrix ud1, matrix ud2)
{
    matrix y;

    // Obliczamy argument funkcji testowej.
    //
    // arg = pi * sqrt((x1/pi)^2 + (x2/pi)^2)
    //
    // Matematycznie jest to po prostu:
    //
    // arg = sqrt(x1^2 + x2^2)
    //
    // czyli odległość punktu od początku układu.

    double arg = M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2));

    // Oryginalna funkcja celu:
    //
    // f(x1,x2) = sin(arg) / arg

    y = sin(arg) / arg;

	// PIERWSZE OGRANICZENIE
//
// g1(x1) = -x1 + 1 <= 0
//
// Jest ono równoważne warunkowi x1 >= 1.
//
// Jeżeli x1 < 1, dodajemy karę.

if (-x(0) + 1 > 0)
    y = y + ud2 * pow(-x(0) + 1, 2);

// DRUGIE OGRANICZENIE
//
// g2(x2) = -x2 + 1 <= 0
//
// Jest ono równoważne warunkowi x2 >= 1.

if (-x(1) + 1 > 0)
    y = y + ud2 * pow(-x(1) + 1, 2);

// TRZECIE OGRANICZENIE
//
// g3(x1,x2) = sqrt(x1^2+x2^2) - a <= 0
//
// Oznacza to, że punkt musi znajdować się
// wewnątrz koła o promieniu a.
//
// ud1 - wartość parametru a
// ud2 - aktualny współczynnik kary c

if (norm(x) - ud1 > 0)
    y = y + ud2 * pow(norm(x) - ud1, 2);

// Zwracamy funkcję celu powiększoną
// o kary za naruszenie ograniczeń.

return y;
}

matrix ff3Tb(matrix x, matrix ud1, matrix ud2)
{
    matrix y;

    // Argument funkcji testowej.
    //
    // Uwaga:
    // można 3.14 zamiast M_PI.

    double arg = 3.14 * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2));

    // Oryginalna funkcja celu.

    y = sin(arg) / arg;

    // PIERWSZE OGRANICZENIE
    //
    // x1 >= 1
    //
    // Jeżeli punkt nie spełnia ograniczenia,
	// przypisujemy bardzo dużą wartość funkcji celu.

	if (-x(0) + 1 > 0)
		y = 1e10;
	else
		// Wewnątrz obszaru dopuszczalnego
		// dodajemy wartość funkcji bariery.
		//
		// Ponieważ g1(x) < 0:
		// -ud2/g1(x) > 0

		y = y - ud2 / (-x(0) + 1);


	// DRUGIE OGRANICZENIE
	//
	// x2 >= 1

	if (-x(1) + 1 > 0)
		y = 1e10;
	else
		y = y - ud2 / (-x(1) + 1);

	// TRZECIE OGRANICZENIE
	//
	// sqrt(x1^2+x2^2) <= a

	if (norm(x) - ud1 > 0)
		y = 1e10;
	else
		y = y - ud2 / (norm(x) - ud1);

	return y;

}

matrix ff3R(matrix x, matrix ud1, matrix ud2)
{
    matrix y;

    // x(0) - początkowa prędkość pozioma v0x
    // x(1) - prędkość kątowa omega
    //
    // Warunki początkowe:
    //
    // x(0) = 0 m
    // vx(0) = v0x
    // y(0) = 100 m
    // vy(0) = 0 m/s

    matrix Y0(4, new double[4] { 0, x(0), 100, 0 });

    // Rozwiązujemy równania ruchu piłki.
    //
    // Czas symulacji: 0-7 s
    // Krok czasowy: 0.01 s
    //
	// Prędkość kątową przekazujemy przez x(1).

	matrix* Y = solve_ode(df3, 0, 0.01, 7, Y0, ud1, x(1));

	int n = get_len(Y[0]);

	// i50 - indeks punktu, dla którego y jest najbliższe 50 m
	// i0 - indeks punktu, dla którego y jest najbliższe 0 m

	int i50 = 0, i0 = 0;

	for (int i = 0; i < n; ++i)
	{
		// Szukamy położenia piłki na wysokości około 50 m.

		if (abs(Y[1](i, 2) - 50) < abs(Y[1](i50, 2) - 50))
			i50 = i;

		// Szukamy punktu uderzenia w podłoże.

		if (abs(Y[1](i, 2)) < abs(Y[1](i0, 2)))
			i0 = i;
	}

	// Funkcja celu:
	//
	// Minimalizujemy ujemny zasięg.
	// xend to położenie poziome przy y ≈ 0.

	y = -Y[1](i0, 0);

	// OGRANICZENIE 1
	//
	// v0x należy do przedziału [-10,10].
	//
	// Jeżeli wartość bezwzględna v0x przekroczy 10,
	// dodajemy karę zewnętrzną.

	if (abs(x(0)) - 10 > 0)
		y = y + ud2 * pow(abs(x(0)) - 10, 2);

	// OGRANICZENIE 2
	//
	// omega należy do przedziału [-10,10].

	if (abs(x(1)) - 10 > 0)
		y = y + ud2 * pow(abs(x(1)) - 10, 2);

	// OGRANICZENIE 3
	//
	// Na wysokości y = 50 m piłka powinna znajdować się
	// w odległości nie większej niż 2 m od punktu x = 5 m.
	//
	// Zatem:
	//
	// |x50 - 5| <= 2
	//
	// czyli:
	//
	// 3 <= x50 <= 7
	//
	// Jeżeli warunek jest niespełniony,
	// dodajemy karę zewnętrzną.

	if (abs(Y[1](i50, 0) - 5) - 2 > 0)
		y = y + ud2 * pow(abs(Y[1](i50, 0) - 5) - 2, 2);

	// Zwracamy ujemny zasięg powiększony o ewentualne kary.

	return y;
}
matrix df3(double t, matrix Y, matrix ud1, matrix ud2)
{
    // Parametry fizyczne modelu.
    //
    // C - współczynnik oporu powietrza
    // r - promień piłki [m]
    // m - masa piłki [kg]
    // ro - gęstość powietrza [kg/m^3]
    // g - przyspieszenie ziemskie [m/s^2]

    double C = 0.47, r = 0.12, m = 0.6, ro = 1.2, g = 9.81;

    // Pole powierzchni przekroju piłki.
    //
    // S = pi * r^2

    double S = M_PI * r * r,

		// Siła oporu w kierunku poziomym.
		//
		// Dx = 0.5 * C * ro * S * vx * |vx|

		Dx = 0.5 * C * ro * S * Y(1) * abs(Y(1)),

		// Siła oporu w kierunku pionowym.
		//
		// Dy = 0.5 * C * ro * S * vy * |vy|

		Dy = 0.5 * C * ro * S * Y(3) * abs(Y(3)),

		// Siła Magnusa w kierunku poziomym.
		//
		// ud2 - prędkość kątowa omega
		//
		// FMx = pi * ro * vy * omega * r^3

		FMx = M_PI * ro * Y(3) * m2d(ud2) * pow(r, 3),

		// Siła Magnusa w kierunku pionowym.
		//
		// FMy = pi * ro * vx * omega * r^3

		FMy = M_PI * ro * Y(1) * m2d(ud2) * pow(r, 3);

    // Tworzymy wektor czterech pochodnych.

    matrix dY(4, 1);

    // POCHODNA POŁOŻENIA POZIOMEGO
    //
    // dx/dt = vx

    dY(0) = Y(1);

    // POCHODNA PRĘDKOŚCI POZIOMEJ
    //
    // m * dvx/dt = -Dx - FMx

    dY(1) = (-Dx - FMx) / m;

    // POCHODNA POŁOŻENIA PIONOWEGO
    //
    // dy/dt = vy

    dY(2) = Y(3);

    // POCHODNA PRĘDKOŚCI PIONOWEJ
    //
    // m * dvy/dt = -m*g - Dy - FMy

    dY(3) = (-m * g - Dy - FMy) / m;

    return dY;
}








