// nao utilizei a forma normalizada dos vetores
// falta a classe do cubo

import java.lang.Math;
class Main {
  public static void main(String[] args) {
    Ponto C = new Ponto(0, 0, 0);
    Vetor n = new Vetor(0, 1, 0);
    double H = 4;
    double R = 4;
    Ponto P0 = new Ponto(3, 0, 3);
    Vetor d = new Vetor(-1, 0.6 , -1);
    Reta X = new Reta(P0, d);
    Cone Y = new Cone(C, R, H, n);
    PontosInt P = Y.IntConReta(X);
    P.PrintPontos();
  }
}


// --------- Ponto ----------
class Ponto {
  double X, Y, Z;

  Ponto(double X, double Y, double Z){
    this.X = X;
    this.Y = Y;
    this.Z = Z;
  }

  void setPonto(double X1, double Y1, double Z1){
    this.X = X1;
    this.Y = Y1;
    this.Z = Z1;
  }

  Ponto SomaPonto(Ponto p){ // soma entre dois pontos
    double XSoma = this.X + p.X;
    double YSoma = this.Y + p.Y;
    double ZSoma = this.Z + p.Z;
    Ponto PSoma = new Ponto(XSoma, YSoma, ZSoma);
    return PSoma;
  }

  Vetor SubPonto(Ponto p){ // subtracao entre dois pontos, resulta num vetor
    double XSub = this.X - p.X;
    double YSub = this.Y - p.Y;
    double ZSub = this.Z - p.Z;
    Vetor PSub = new Vetor(XSub, YSub, ZSub);
    return PSub;
  }

  double Distancia(Ponto p){
    double difX = Math.abs(this.X - p.X);
    double difY = Math.abs(this.Y - p.Y);
    double difZ = Math.abs(this.Z - p.Z);
    double DistXY = Math.sqrt((difX*difX) + (difY*difY));
    double DistXYZ = Math.sqrt((DistXY*DistXY) + (difZ*difZ));
    return DistXYZ;
  }

  public String toString(){
    return "(" + this.X + ", " + this.Y + ", " + this.Z + ")";
  }
}


// -------- Pontos de Intersecao ----------- (classe que eu criei pq nao consegui fazer nada melhor)
class PontosInt{ // Classe feita para guardar os pontos de intersecao
  Ponto pontos[];
  int N;

  PontosInt(){
    this.pontos = new Ponto[2];
    this.N = 0;
  }

  PontosInt(Ponto p){
    this.pontos = new Ponto[2];
    this.pontos[0] = p;
    this.N = 1;
  }

  PontosInt(Ponto p1, Ponto p2){
    this.pontos = new Ponto[2];
    this.pontos[0] = p1;
    this.pontos[1] = p2;
    this.N = 2;
  }

  void addPonto(Ponto p){
    this.pontos[this.N] = p;
    this.N += 1;
  }

  public void PrintPontos(){
    int i = 0;
    for (i = 0; i < N; i++){
      System.out.println(pontos[i].toString());
    }
    if (i == 0){
      System.out.println("Nao ha intersecao");
    }
  }
}


// ---------- Reta ----------
class Reta{
  Ponto p0; // Ponto inicial da reta
  double t; // Escalar
  Vetor d; // Vetor unitario de sentido e direcao da reta

  Reta(Ponto p0, Vetor d){
    this.p0 = p0;
    this.d = d;
  }

  Ponto PnaReta(double t){ // P(t) = p0 + t*d
    Vetor V = this.d.ProdEscalar(t);
    Ponto P =  V.SomaPtVt(this.p0);
    return P;
  }
}


// --------- Vetor ----------
class Vetor{
  double X, Y, Z;

  Vetor(double X, double Y, double Z){
    this.X = X;
    this.Y = Y;
    this.Z = Z;
  }

  Vetor ProdEscalar(double escalar){ // produto entre um vetor e um escalar
    double X1 = this.X * escalar;
    double Y1 = this.Y * escalar;
    double Z1 = this.Z * escalar;
    Vetor Prod = new Vetor(X1, Y1, Z1);
    return Prod;
  }

  double ProdVetorial(Vetor v){ // produto entre dois vetores
    double Prod = (this.X * v.X) + (this.Y * v.Y) + (this.Z * v.Z);
    return Prod;
  }

  double ProdVetPt(Ponto p){ // produto entre vetor e ponto
    double Prod = (this.X * p.X) + (this.Y * p.Y) + (this.Z * p.Z);
    return Prod;
  }

  Ponto SomaPtVt(Ponto p){ // soma entre um ponto e um vetor
    double X1 = this.X + p.X;
    double Y1 = this.Y + p.Y;
    double Z1 = this.Z + p.Z;
    return new Ponto(X1, Y1, Z1);
  }

  Ponto SubPtVt(Ponto p){ // subtracao entre um ponto e um vetor
    double X1 = p.X - this.X;
    double Y1 = p.Y - this.Y;
    double Z1 = p.Z - this.Z;
    return new Ponto(X1, Y1, Z1);
  }

  Vetor SubVetor(Vetor v){ // subtracao entre dois vetores
    double XSub = this.X - v.X;
    double YSub = this.Y - v.Y;
    double ZSub = this.Z - v.Z;
    Vetor Sub = new Vetor(XSub, YSub, ZSub);
    return Sub;
  }
}


// ---------- Plano ----------
class Plano{
  Ponto Ppl; // Ponto qualquer do plano
  Vetor n; // Vetor unitario perpendicular ao plano

  Plano(Ponto P0, Vetor n){
    this.Ppl = P0;
    this.n = n;
  }

  Ponto InterRetaPlano(Reta R){ // (P - Ppl) * n = 0
    Vetor PplP0 = this.Ppl.SubPonto(R.p0);
    double T =  PplP0.ProdVetorial(this.n) / this.n.ProdVetorial(R.d);
    Ponto PInt = R.PnaReta(T);
    return PInt;
  }
}


// --------- Muro ---------
class Painel extends Plano{ // Vetor do muro precisa ser paralelo ao eixo Z (ou outro eixo, nesse caso precisa alterar o codigo)
  double L; // L = lado
  Ponto buracos[][];
  int H, V;

  Painel(double L, Vetor n, Ponto C, int h, int v){
    super(new Ponto((C.X - L/2), (C.Y - L/2), C.Z), n); // o ponto de referencia eh o do canto inferior esquerdo e nao o centro do painel
    this.L = L;
    this.buracos = new Ponto[h][v];
    this.H = h;
    this.V = v;
  }
  
  Ponto PnoPainel(int h, int v){ // h = posicao na matriz da esquerda para direita, v = posicao na matriz de baixo para cima, L/2H+(h-1)*L/H
    Ponto p = new Ponto(Ppl.X + (this.L/2*this.H) + (h-1)*L/H, Ppl.Y + (this.L/2*this.V) + (h-1)*L/H, Ppl.Z);
    return p;
  }
}


// -------- Esfera ---------
class Esfera{
  Ponto C; // Ponto do centro da esfera
  double R; // R = raio

  Esfera(Ponto C, double R){
    this.C = C;
    this.R = R;
  }

  double CalculoDeltaEsf(Reta R){
    double a = R.d.ProdVetorial(R.d); // a = (d*d)
    double b = (R.p0.SubPonto(this.C)).ProdVetorial(R.d); // b = ((p0-C)*d)
    double c = (R.p0.SubPonto(this.C)).ProdVetorial(R.p0.SubPonto(this.C)) - (this.R * this.R); // c = ((p0-C)*(p0-C)-(R**2))
    double Delta = (b*b) - (a*c);
    return Delta;
  }

  PontosInt IntEsfReta(Reta R){ // (P-C)*(P-C) = R**2
    double a = R.d.ProdVetorial(R.d);
    double b = (R.p0.SubPonto(this.C)).ProdVetorial(R.d);
    double Delta = this.CalculoDeltaEsf(R);
    if (Delta < 0) {
      return null;
    }
    if (Delta == 0){
      double T = -b/a;
      Ponto p1 = R.PnaReta(T);
      return new PontosInt(p1);
    }
    if (Delta > 0){
      double T = (-b + Math.sqrt(Delta)) / a;
      Ponto p1 = R.PnaReta(T);
      T = (-b - Math.sqrt(Delta)) / a;
      Ponto p2 = R.PnaReta(T);
      return new PontosInt(p1, p2);
    }
    return null;
  }
}


//  --------- Cilindro ---------
class Cilindro{
  Ponto B; // Ponto da base do cilindro
  double R, H; // R = raio; H = altura
  Vetor u; // Vetor unitario da direcao e sentido do cilindro

  Cilindro(Ponto B, Vetor u, double R, double H){
    this.B = B;
    this.u = u;
    this.R = R;
    this.H = H;
  }

  double CalculoDeltaCil(Reta r){
    Vetor w = r.d.SubVetor(this.u.ProdEscalar(r.d.ProdVetorial(this.u))); // w = d - (d*u)*u
    Vetor v = r.p0.SubPonto(this.B).SubVetor(this.u.ProdEscalar(r.p0.SubPonto(this.B).ProdVetorial(this.u))); // v = (p0-B)-((p0-B)*u)*u
    double a = w.ProdVetorial(w); // a = (w*w)
    double b = v.ProdVetorial(w); // b = (v*w)
    double c = v.ProdVetorial(v) - (this.R * this.R); // c = (v*v) - R**2
    double Delta = (b*b) - (a*c);
    return Delta;
  }

  PontosInt IntCilReta(Reta r){ 
    Vetor w = r.d.SubVetor(this.u.ProdEscalar(r.d.ProdVetorial(this.u)));
    Vetor v = r.p0.SubPonto(this.B).SubVetor(this.u.ProdEscalar(r.p0.SubPonto(this.B).ProdVetorial(this.u)));
    double a = w.ProdVetorial(w);
    double b = v.ProdVetorial(w);
    double Delta = this.CalculoDeltaCil(r);
    if (Delta < 0){
      if (a == 0){
        Plano pl = new Plano(this.B, this.u);
        Ponto IntPlanoReta = pl.InterRetaPlano(r);
        double Dist = IntPlanoReta.Distancia(this.B);
        PontosInt PtsValidos = new PontosInt();
        if (Dist < this.R){
          Ponto topo = r.d.ProdEscalar(H).SomaPtVt(IntPlanoReta);
          PtsValidos.addPonto(IntPlanoReta);
          PtsValidos.addPonto(topo);
        }
        return PtsValidos;
      }
      return null;
    }
    if (Delta == 0){
      double T = -b / a;
      Ponto p1 = r.PnaReta(T);
      double Verif = p1.SubPonto(this.B).ProdVetorial(this.u);
      if (Verif >= 0 && Verif <= H){
        return new PontosInt(p1);
      }
    return null;
    }
    if (Delta > 0){
      boolean p1V = false, p2V = false;
      double T = (-b + Math.sqrt(Delta)) / a;
      Ponto p1 = r.PnaReta(T);
      double Verif1 = p1.SubPonto(this.B).ProdVetorial(this.u);
      PontosInt PtsValidos = new PontosInt();
      if (Verif1 >= 0 && Verif1 <= H){
        PtsValidos.addPonto(p1);
        p1V = true;
      }
      T = (-b - Math.sqrt(Delta)) / a;
      Ponto p2 = r.PnaReta(T);
      double Verif2 = p2.SubPonto(this.B).ProdVetorial(this.u);
      if (Verif2 >= 0 && Verif2 <= H){
        PtsValidos.addPonto(p2);
        p2V = true;
      }
      if ((p1V && !(p2V)) || !((p1V) && p2V)){ // Apenas um pt valido
        Plano pl = new Plano(this.B, this.u);
        Ponto IntPlReta = pl.InterRetaPlano(r);
        double Dist = IntPlReta.Distancia(B);
        if (Dist < this.R){
          PtsValidos.addPonto(IntPlReta);
        }
        else{
          Ponto Int = r.d.ProdEscalar(H).SomaPtVt(IntPlReta);
          PtsValidos.addPonto(Int);
        }
      }
      if (!p1V && !p2V && ((Verif1 < 0 && Verif2 > H)||(Verif1 > H && Verif2 < 0))){ // Dois invalidos e um acima e outro abaixo do cilindro
        Plano pl = new Plano(this.B, this.u);
        Ponto IntPlaReta = pl.InterRetaPlano(r);
        Ponto IntRetaTopo = r.d.ProdEscalar(H).SomaPtVt(IntPlaReta);
        PtsValidos.addPonto(IntPlaReta);
        PtsValidos.addPonto(IntRetaTopo);
      }
      return PtsValidos;
    }
    return null;
  }
}


// -------- Cone ---------
class Cone{ // cos**2(X) = H**2 / H**2 + R**2
  Vetor n; // Vetor direcional do cone
  Ponto C, V; // C = centro da base, V = vertice do cone
  double H, R; // H = altura, R = raio

  Cone(Ponto C, Ponto V, double R, double H, Vetor n){
    this.n = n;
    this.H = H;
    this.R = R;
    this.C = C;
    this.V = V;
  }

  Cone(Ponto C, double R, double H, Vetor n){
    this.C = C;
    this.n = n;
    this.R = R;
    this.H = H;
    this.V = n.ProdEscalar(H).SomaPtVt(C);
  }

  Cone(double R, double H, Vetor n, Ponto V){
    this.V = V;
    this.H = H;
    this.R = R;
    this.n = n;
    this.C = n.ProdEscalar(H).SubPtVt(V);
  }

  double CalculoDeltaCon(Reta r){
    Vetor v = V.SubPonto(r.p0);
    double a = (r.d.ProdVetorial(this.n))*(r.d.ProdVetorial(this.n)) - (r.d.ProdVetorial(r.d) * (H*H/((H*H)+(R*R)))); // (d*n)**2-(d*d)*cos**2(X)
    double b = v.ProdVetorial(r.d)*(H*H/((H*H)+(R*R)))-(v.ProdVetorial(n) * r.d.ProdVetorial(n)); // (v*d)*cos**2(X) - (v*n)*(d*n)
    double c = v.ProdVetorial(n)*v.ProdVetorial(n)-v.ProdVetorial(v)*(H*H/((H*H)+(R*R))); // (v*n)**2 - (v*v)*cos**2(X)
    double Delta = (b*b) - (a*c);
    return Delta;
  }

  PontosInt IntConReta(Reta r){
    Vetor v = V.SubPonto(r.p0);
    double a = (r.d.ProdVetorial(this.n))*(r.d.ProdVetorial(this.n)) - (r.d.ProdVetorial(r.d) * (H*H/((H*H)+(R*R))));
    double b = v.ProdVetorial(r.d)*(H*H/((H*H)+(R*R)))-(v.ProdVetorial(n) * r.d.ProdVetorial(n));
    double c = v.ProdVetorial(n)*v.ProdVetorial(n)-v.ProdVetorial(v)*(H*H/((H*H)+(R*R)));
    double Delta = this.CalculoDeltaCon(r);
    if (a == 0){
      double T = -(c / (2 * b));
      Ponto p1 = r.PnaReta(T);
      if (V.SubPonto(p1).ProdVetorial(n) >= 0 && V.SubPonto(p1).ProdVetorial(n) <= H){
        return new PontosInt(p1);
      }
      return null;
    }
    if (Delta < 0){
      return null;
    }
    if (Delta == 0){
      double T = -b / a;
      Ponto p1 = r.PnaReta(T);
      if (V.SubPonto(p1).ProdVetorial(n) >= 0 && V.SubPonto(p1).ProdVetorial(n) <= H){
        return new PontosInt(p1);
      }
      return null;
    }
    if (Delta > 0){
      PontosInt IntRetaCon = new PontosInt();
      boolean PV1 = false, PV2 = false;
      double T1 = (-b + Math.sqrt(Delta)) / a;
      Ponto p1 = r.PnaReta(T1);
      if (V.SubPonto(p1).ProdVetorial(n) >= 0 && V.SubPonto(p1).ProdVetorial(n) <= H){
        PV1 = true;
        IntRetaCon.addPonto(p1);
      }
      double T2 = (-b - Math.sqrt(Delta)) / a;
      Ponto p2 = r.PnaReta(T2);
      if (V.SubPonto(p2).ProdVetorial(n) >= 0 && V.SubPonto(p2).ProdVetorial(n) <= H){
        PV2 = true;
        IntRetaCon.addPonto(p2);
      }
      Plano pl = new Plano(this.C, this.n);
      Ponto IntRetaPl = pl.InterRetaPlano(r);
      if(IntRetaPl.Distancia(this.C) < this.R){
        IntRetaCon.addPonto(IntRetaPl);
      }
      return IntRetaCon;
    }
    return null;
  }
}