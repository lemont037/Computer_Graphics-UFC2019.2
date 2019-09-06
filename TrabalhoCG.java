import java.lang.Math;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.text.DecimalFormat;

class Main {
  public static void main(String[] args) {
    Espaco E = new Espaco();
    E.addObjeto(new Cilindro(new Ponto(0,-2,-10), new Vetor(0,1,0), 0.5, 2, new RGB(139,69,19), 1));
    E.addObjeto(new Cone(new Ponto(0,0,-10), 3, 8, new Vetor(0,1,0), new RGB(0,255,127), 2));
    E.addObjeto(new Cubo(new Ponto(0, 1,-20), 6, new RGB(178,34,34), 3));
    E.addObjeto(new Cubo(new Ponto(0, 7,-20), 6, new RGB(178,34,34), 4));
    E.addObjeto(new Cubo(new Ponto(0, 13,-20), 6, new RGB(178,34,34), 5));
    //Ponto LookAt = new Ponto();
    //Ponto ViewUp = new Ponto();
    Ponto observador = new Ponto(0,0,0);
    //Vetor K = observador.SubPonto(LookAt);
    //Vetor k = K.ProdEscalar(1/sqrt(K.ProdEscalar(K));
    Muro painel = new Muro(4, new Vetor(0,0,1),new Ponto(0,0,4), 400, 400);
    Imagem I = new Imagem(400,400);
    for (int i = 0; i < 400; i++){
      for (int j = 0; j < 400; j++){
        Ponto Aux = painel.PnoMuro(i, j);
        Vetor n = Aux.SubPonto(observador);
        Reta R = new Reta(observador, n);
        PontosInt Pfinal = new PontosInt();
        for (int k = 0; k < E.N; k++){
          PontosInt p = E.objetos[k].InterReta(R);
          Pfinal = Pfinal.Uniao(p);
        }
        if (Pfinal.N != 0){
          FigurasGeo prim = E.getFigPorId(Pfinal.getPonto(0).id);
          Pfinal.setPrimObj(prim);
          Pixel Pix = new Pixel(Pfinal);
          I.addPixel(Pix, i, j);

        }
        else{
          I.pixelNulo(i, j);
        }
      }
    }
    System.out.println(I.pixels[200][150].printOrdem());
    System.out.println(I.pixels[200][150].printPontos());
    Quadro quadro = new Quadro(I);
    //Janela janela = new Janela(E, painel, I);
  }
}

//---------- Quadro --------------
class Quadro extends Frame{
  Quadro(Imagem I){
    add(new DrawingComponent(I));
    setVisible(true);
    setSize(I.H,I.V);
    setLocation(150,150);
    addWindowListener(new WindowAdapter(){
      public void windowClosing(WindowEvent f){
        dispose();  }    });
  }
}
  


// --------- Ponto ----------
class Ponto {
  double X, Y, Z, T;
  int id; // figura ao qual o ponto pertence
  RGB cor;

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

  void setT(double T){
    this.T = T;
  }

  void setID(int id){
    this.id = id;
  }

  void setCor(RGB cor){
    this.cor = cor;
  }

  double getCoordX(){
    return this.X;
  } 

  double getCoordY(){
    return this.Y;
  }

  double getCoordZ(){
    return this.Z;
  }

  Ponto SomaPonto(Ponto p){
    double XSoma = this.X + p.X;
    double YSoma = this.Y + p.Y;
    double ZSoma = this.Z + p.Z;
    Ponto PSoma = new Ponto(XSoma, YSoma, ZSoma);
    return PSoma;
  }

  Vetor SubPonto(Ponto p){
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
    DecimalFormat df = new DecimalFormat("0.000");
    return "(" + df.format(this.X) + ", " + df.format(this.Y) + ", " + df.format(this.Z) + ")";
  }

  String Id(){
    return Integer.toString(id);
  }
}

// -------- Pontos de Intersecao ----------- (classe que eu criei pq nao consegui fazer nada melhor)
class PontosInt{ // Classe feita para guardar os pontos de intersecao
  Ponto pontos[];
  int N;
  FigurasGeo PrimObj;

  PontosInt(){
    this.pontos = new Ponto[50];
    this.N = 0;
  }

  PontosInt(Ponto p){
    this.pontos = new Ponto[50];
    this.pontos[0] = p;
    this.N = 1;
  }

  PontosInt(Ponto p1, Ponto p2){
    this.pontos = new Ponto[50];
    this.pontos[0] = p1;
    this.pontos[1] = p2;
    this.N = 2;
  }

  Ponto getPonto(int posicao){
    if (this.N > 0){
      return this.pontos[posicao];
    }
    return null;
  }

  void setPrimObj(FigurasGeo p){
    this.PrimObj = p;
  }

  void addPonto(Ponto p){
    if (this.N != 0 && pontos[0].T < p.T){
      for (int i = 0; i < N; i++){
        pontos[i+1] = pontos[i];
      }
      this.pontos[0] = p;
    }
    else{
      this.pontos[this.N] = p;
    }
    this.N += 1;
  }

  PontosInt Uniao(PontosInt p){
    if (this.N == 0){
      return p;
    }
    if (p.N == 0){
      return this;
    }
    PontosInt uniao = new PontosInt();
    for (int i = 0; i < p.N; i++){
      uniao.addPonto(p.getPonto(i));
    }
    for (int i = 0; i < this.N; i++){
      uniao.addPonto(this.getPonto(i));
    }
    return uniao;
  }

  String ObjetosId(){
    if (this.N == 0){
      return "Nao ha intersecao";
    }
    String Ids = new String(pontos[0].Id());
    for (int i = 1; i < this.N; i++){
      String aux1 = pontos[i-1].Id();
      String aux2 = pontos[i].Id();
      if (!(aux1.equals(aux2))){
        Ids = Ids.concat(" " + pontos[i].Id());
      }
    }
    return Ids;
  }

  public String toString(){
    if (this.N == 0){
      return "Nao ha intersecao";
    }
    String pts = new String();
    for (int i = 0; i< this.N; i++){
      pts = pts.concat(pontos[i].toString() + "\n");
    }
    return pts;
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

  Vetor ProdEscalar(double escalar){
    double X1 = this.X * escalar;
    double Y1 = this.Y * escalar;
    double Z1 = this.Z * escalar;
    Vetor Prod = new Vetor(X1, Y1, Z1);
    return Prod;
  }

  double ProdVetorial(Vetor v){
    double Prod = (this.X * v.X) + (this.Y * v.Y) + (this.Z * v.Z);
    return Prod;
  }

  double ProdVetPt(Ponto p){
    double Prod = (this.X * p.X) + (this.Y * p.Y) + (this.Z * p.Z);
    return Prod;
  }

  Ponto SomaPtVt(Ponto p){
    double X1 = this.X + p.X;
    double Y1 = this.Y + p.Y;
    double Z1 = this.Z + p.Z;
    return new Ponto(X1, Y1, Z1);
  }

  Ponto SubPtVt(Ponto p){
    double X1 = p.X - this.X;
    double Y1 = p.Y - this.Y;
    double Z1 = p.Z - this.Z;
    return new Ponto(X1, Y1, Z1);
  }

  Vetor SubVetor(Vetor v){
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

  Plano(Ponto P, Vetor n){
    this.Ppl = P;
    this.n = n;
  }

  Ponto InterRetaPlano(Reta R){ // (P - Ppl) * n = 0
    if (this.n.ProdVetorial(R.d) == 0){
      return null;
    }
    Vetor PplP0 = this.Ppl.SubPonto(R.p0);
    double T =  PplP0.ProdVetorial(this.n) / this.n.ProdVetorial(R.d);
    Ponto PInt = R.PnaReta(T);
    PInt.setT(T);
    return PInt;
  }
}

// --------- Muro ---------
class Muro extends Plano{ // Vetor do muro precisa ser paralelo ao eixo Z
  double L; // L = lado
  Ponto buracos[][];
  int H, V;

  Muro(double L, Vetor n, Ponto C, int h, int v){
    super(new Ponto((C.X - L/2), (C.Y + L/2), C.Z), n);
    this.L = L;
    this.buracos = new Ponto[h][v];
    this.H = h;
    this.V = v;
  }

  Ponto PnoMuro(int h, int v){ // Calculo das coords X e Y do ponto no muro, L/2H+(h-1)*L/H
    Ponto p = new Ponto(this.Ppl.X + (this.L/(2*this.H)) + (h)*this.L/this.H, this.Ppl.Y - (this.L/(2*this.V)) - (v)*this.L/this.V, this.Ppl.Z);
    this.buracos[h][v] = p;
    return p;
  }
}

// -------- Esfera ---------
class Esfera extends FigurasGeo{
  Ponto C; // Ponto do centro da esfera
  double R; // R = raio

  Esfera(Ponto C, double R, RGB cor, int id){
    this.C = C;
    this.R = R;
    this.cor = cor;
    this.id = id;
    this.acertado = false;
  }

  double CalculoDeltaEsf(Reta R){
    double a = R.d.ProdVetorial(R.d); // a = (d*d)
    double b = (R.p0.SubPonto(this.C)).ProdVetorial(R.d); // b = ((p0-C)*d)
    double c = (R.p0.SubPonto(this.C)).ProdVetorial(R.p0.SubPonto(this.C)) - (this.R * this.R); // c = ((p0-C)*(p0-C)-(R**2))
    double Delta = (b*b) - (a*c);
    return Delta;
  }

  PontosInt InterReta(Reta r){ // (P-C)*(P-C) = R**2
    double a = r.d.ProdVetorial(r.d);
    double b = (r.p0.SubPonto(this.C)).ProdVetorial(r.d);
    double Delta = this.CalculoDeltaEsf(r);
    if (Delta < 0) {
      new PontosInt();
    }
    if (Delta == 0){
      double T = -b/a;
      Ponto p1 = r.PnaReta(T);
      p1.setT(T);
      p1.setCor(this.cor);
      p1.setID(id);
      return new PontosInt(p1);
    }
    if (Delta > 0){
      double T = (-b + Math.sqrt(Delta)) / a;
      Ponto p1 = r.PnaReta(T);
      p1.setT(T);
      p1.setID(this.id);
      p1.setCor(this.cor);
      T = (-b - Math.sqrt(Delta)) / a;
      Ponto p2 = r.PnaReta(T);
      p2.setT(T);
      p2.setID(this.id);
      p2.setCor(this.cor);
      return new PontosInt(p1, p2);
    }
    return new PontosInt();
  }
}


//  --------- Cilindro ---------
class Cilindro extends FigurasGeo{
  Ponto B; // Ponto da base do cilindro
  double R, H; // R = raio; H = altura
  Vetor u; // Vetor unitario da direcao e sentido do cilindro

  Cilindro(Ponto B, Vetor u, double R, double H, RGB cor, int id){
    this.B = B;
    this.u = u;
    this.R = R;
    this.H = H;
    this.cor = cor;
    this.id = id;
    this.acertado = false;
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

  PontosInt InterReta(Reta r){ 
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
          topo.setT(H + IntPlanoReta.T);
          topo.setCor(cor);
          topo.setID(id);
          IntPlanoReta.setID(id);
          IntPlanoReta.setCor(cor);
          PtsValidos.addPonto(IntPlanoReta);
          PtsValidos.addPonto(topo);
        }
        return PtsValidos;
      }
      return new PontosInt();
    }
    if (Delta == 0){
      double T = -b / a;
      Ponto p1 = r.PnaReta(T);
      double Verif = p1.SubPonto(this.B).ProdVetorial(this.u);
      if (Verif >= 0 && Verif <= H){
        p1.setT(T);
        p1.setID(this.id);
        p1.setCor(cor);
        return new PontosInt(p1);
      }
    return new PontosInt();
    }
    if (Delta > 0){
      boolean p1V = false, p2V = false;
      double T = (-b + Math.sqrt(Delta)) / a;
      Ponto p1 = r.PnaReta(T);
      double Verif1 = p1.SubPonto(this.B).ProdVetorial(this.u);
      PontosInt PtsValidos = new PontosInt();
      if (Verif1 >= 0 && Verif1 <= H){
        p1.setT(T);
        p1.setID(id);
        p1.setCor(cor);
        PtsValidos.addPonto(p1);
        p1V = true;
      }
      T = (-b - Math.sqrt(Delta)) / a;
      Ponto p2 = r.PnaReta(T);
      double Verif2 = p2.SubPonto(this.B).ProdVetorial(this.u);
      if (Verif2 >= 0 && Verif2 <= H){
        p2.setT(T);
        p2.setID(id);
        p2.setCor(cor);
        PtsValidos.addPonto(p2);
        p2V = true;
      }
      if ((p1V && !(p2V)) || (!(p1V) && p2V)){ // Apenas um pt valido
        Plano pl = new Plano(this.B, this.u);
        Ponto IntPlReta = pl.InterRetaPlano(r);
        double Dist = IntPlReta.Distancia(B);
        if (Dist < this.R){
          IntPlReta.setID(id);
          IntPlReta.setCor(cor);
          PtsValidos.addPonto(IntPlReta);
        }
        else{
          Ponto Int = r.d.ProdEscalar(H).SomaPtVt(IntPlReta);
          Int.setT(H + IntPlReta.T);
          Int.setID(id);
          Int.setCor(cor);
          PtsValidos.addPonto(Int);
        }
      }
      if (!(p1V) && !(p2V) && ((Verif1 < 0 && Verif2 > H)||(Verif1 > H && Verif2 < 0))){ // Dois invalidos e um acima e outro abaixo do cilindro
        Plano pl = new Plano(this.B, this.u);
        Ponto IntPlaReta = pl.InterRetaPlano(r);
        Ponto IntRetaTopo = r.d.ProdEscalar(H).SomaPtVt(IntPlaReta);
        IntRetaTopo.setT(H + IntPlaReta.T);
        IntRetaTopo.setID(id);
        IntRetaTopo.setCor(cor);
        IntPlaReta.setID(id);
        IntPlaReta.setCor(cor);
        PtsValidos.addPonto(IntPlaReta);
        PtsValidos.addPonto(IntRetaTopo);
      }
      return PtsValidos;
    }
    return new PontosInt();
  }
}


// -------- Cone ---------
class Cone extends FigurasGeo{ // cos**2(X) = H**2 / H**2 + R**2
  Vetor n; // Vetor direcional do cone
  Ponto C, V; // C = centro da base, V = vertice do cone
  double H, R; // H = altura, R = raio

  Cone(Ponto C, Ponto V, double R, double H, Vetor n, RGB cor, int id){
    this.n = n;
    this.H = H;
    this.R = R;
    this.C = C;
    this.V = V;
    this.cor = cor;
    this.id = id;
    this.acertado = false;
  }

  Cone(Ponto C, double R, double H, Vetor n, RGB cor, int id){
    this.C = C;
    this.n = n;
    this.R = R;
    this.H = H;
    this.V = n.ProdEscalar(H).SomaPtVt(C);
    this.cor = cor;
    this.id = id;
    this.acertado = false;
  }

  Cone(double R, double H, Vetor n, Ponto V, RGB cor, int id){
    this.V = V;
    this.H = H;
    this.R = R;
    this.n = n;
    this.C = n.ProdEscalar(H).SubPtVt(V);
    this.cor = cor;
    this.id = id;
    this.acertado = false;
  }

  double CalculoDeltaCon(Reta r){
    Vetor v = V.SubPonto(r.p0);
    double a = (r.d.ProdVetorial(this.n))*(r.d.ProdVetorial(this.n)) - (r.d.ProdVetorial(r.d) * (H*H/((H*H)+(R*R)))); // (d*n)**2-(d*d)*cos**2(X)
    double b = v.ProdVetorial(r.d)*(H*H/((H*H)+(R*R)))-(v.ProdVetorial(n) * r.d.ProdVetorial(n)); // (v*d)*cos**2(X) - (v*n)*(d*n)
    double c = v.ProdVetorial(n)*v.ProdVetorial(n)-v.ProdVetorial(v)*(H*H/((H*H)+(R*R))); // (v*n)**2 - (v*v)*cos**2(X)
    double Delta = (b*b) - (a*c);
    return Delta;
  }

  PontosInt InterReta(Reta r){
    Vetor v = V.SubPonto(r.p0);
    double a = (r.d.ProdVetorial(this.n))*(r.d.ProdVetorial(this.n)) - (r.d.ProdVetorial(r.d) * (H*H/((H*H)+(R*R))));
    double b = v.ProdVetorial(r.d)*(H*H/((H*H)+(R*R)))-(v.ProdVetorial(n) * r.d.ProdVetorial(n));
    double c = v.ProdVetorial(n)*v.ProdVetorial(n)-v.ProdVetorial(v)*(H*H/((H*H)+(R*R)));
    double Delta = this.CalculoDeltaCon(r);
    if (a == 0){
      double T = -(c / (2 * b));
      Ponto p1 = r.PnaReta(T);
      if (V.SubPonto(p1).ProdVetorial(n) >= 0 && V.SubPonto(p1).ProdVetorial(n) <= H){
        p1.setT(T);
        p1.setCor(cor);
        p1.setID(this.id);
        return new PontosInt(p1);
      }
      return new PontosInt();
    }
    if (Delta < 0){
      return new PontosInt();
    }
    if (Delta == 0){
      double T = -b / a;
      Ponto p1 = r.PnaReta(T);
      if (V.SubPonto(p1).ProdVetorial(n) >= 0 && V.SubPonto(p1).ProdVetorial(n) <= H){
        p1.setT(T);
        p1.setCor(cor);
        p1.setID(this.id);
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
        p1.setT(T1);
        p1.setCor(cor);
        p1.setID(this.id);
        IntRetaCon.addPonto(p1);
      }
      double T2 = (-b - Math.sqrt(Delta)) / a;
      Ponto p2 = r.PnaReta(T2);
      if (V.SubPonto(p2).ProdVetorial(n) >= 0 && V.SubPonto(p2).ProdVetorial(n) <= H){
        PV2 = true;
        p2.setT(T2);
        p2.setCor(cor);
        p2.setID(this.id);
        IntRetaCon.addPonto(p2);
      }
      Plano pl = new Plano(this.C, this.n);
      Ponto IntRetaPl = pl.InterRetaPlano(r);
      if(IntRetaPl.Distancia(this.C) < this.R){
        IntRetaPl.setID(this.id);
        IntRetaPl.setCor(cor);
        IntRetaCon.addPonto(IntRetaPl);
      }
      return IntRetaCon;
    }
    return new PontosInt();
  }
}

// ---------- Cubo ----------
class Cubo extends FigurasGeo{
  double A; // aresta
  Ponto p0, p1, p2, p3, p4, p5, p6, p7, p8;

  Cubo(Ponto p, double A, RGB cor, int id){
    this.p0 = p;
    this.A = A;
    this.cor = cor;
    this.id = id;
    this.acertado = false;
    this.p1 = new Ponto(p.X+A/2, p.Y+A/2, p.Z+A/2);
    this.p2 = new Ponto(p.X+A/2, p.Y+A/2, p.Z-A/2);
    this.p3 = new Ponto(p.X+A/2, p.Y-A/2, p.Z+A/2);
    this.p4 = new Ponto(p.X+A/2, p.Y-A/2, p.Z-A/2);
    this.p5 = new Ponto(p.X-A/2, p.Y+A/2, p.Z+A/2);
    this.p6 = new Ponto(p.X-A/2, p.Y+A/2, p.Z-A/2);
    this.p7 = new Ponto(p.X-A/2, p.Y-A/2, p.Z+A/2);
  }

  PontosInt InterReta(Reta r){
    PontosInt ptsInt = new PontosInt();
    Ponto P1 = new Ponto(p0.X + A/2, p0.Y, p0.Z);
    Ponto P2 = new Ponto(p0.X - A/2, p0.Y, p0.Z);
    Ponto P3 = new Ponto(p0.X, p0.Y + A/2, p0.Z);
    Ponto P4 = new Ponto(p0.X, p0.Y - A/2, p0.Z);
    Ponto P5 = new Ponto(p0.X, p0.Y, p0.Z + A/2);
    Ponto P6 = new Ponto(p0.X, p0.Y, p0.Z - A/2);
    Plano face1 = new Plano(P1, new Vetor(1,0,0));
    Plano face2 = new Plano(P2, new Vetor(1,0,0));
    Plano face3 = new Plano(P3, new Vetor(0,1,0));
    Plano face4 = new Plano(P4, new Vetor(0,1,0));
    Plano face5 = new Plano(P5, new Vetor(0,0,1));
    Plano face6 = new Plano(P6, new Vetor(0,0,1));
    Ponto Int1 = face1.InterRetaPlano(r);
    Ponto Int2 = face2.InterRetaPlano(r);
    Ponto Int3 = face3.InterRetaPlano(r);
    Ponto Int4 = face4.InterRetaPlano(r);
    Ponto Int5 = face5.InterRetaPlano(r);
    Ponto Int6 = face6.InterRetaPlano(r);
    if (Math.abs(Int1.X - P1.X) <= A/2 && Math.abs(Int1.Y - P1.Y) <= A/2 && Math.abs(Int1.Z - P1.Z) <= A/2){
      Int1.setID(id);
      Int1.setCor(this.cor);
      ptsInt.addPonto(Int1);
    }
    if (Math.abs(Int2.X - P2.X) <= A/2 && Math.abs(Int2.Y - P2.Y) <= A/2 && Math.abs(Int2.Z - P2.Z) <= A/2){
      Int2.setID(id);
      Int2.setCor(this.cor);
      ptsInt.addPonto(Int2);
    }
    if (Math.abs(Int3.X - P3.X) <= A/2 && Math.abs(Int3.Y - P3.Y) <= A/2 && Math.abs(Int3.Z - P3.Z) <= A/2){
      Int3.setID(id);
      Int3.setCor(this.cor);
      ptsInt.addPonto(Int3);
    }
    if (Math.abs(Int4.X - P4.X) <= A/2 && Math.abs(Int4.Y - P4.Y) <= A/2 && Math.abs(Int4.Z - P4.Z) <= A/2){
      Int4.setID(id);
      Int4.setCor(this.cor);
      ptsInt.addPonto(Int4);
    }
    if (Math.abs(Int5.X - P5.X) <= A/2 && Math.abs(Int5.Y - P5.Y) <= A/2 && Math.abs(Int5.Z - P5.Z) <= A/2){
      Int5.setID(id);
      Int5.setCor(this.cor);
      ptsInt.addPonto(Int5);
    }
    if (Math.abs(Int6.X - P6.X) <= A/2 && Math.abs(Int6.Y - P6.Y) <= A/2 && Math.abs(Int6.Z - P6.Z) <= A/2){
      Int6.setID(id);
      Int6.setCor(this.cor);
      ptsInt.addPonto(Int6);
    }
    return ptsInt;
  }
}


// ---------- Espaco ----------
class Espaco{
  int N; // N = numero de objetos no Espaco
  FigurasGeo objetos[];

  Espaco(){
    this.N = 0;
    this.objetos = new FigurasGeo[10];
  }

  void addObjeto(FigurasGeo obj){
    objetos[this.N] = obj;
    this.N += 1;
  }

  FigurasGeo getFigPorId(int id){
    for (int i = 0; i < N; i++){
      if (objetos[i].id == id){
        return objetos[i];
      }
    }
    return null;
  }
}


// ---------- Figuras Geometricas ---------
abstract class FigurasGeo{
  RGB cor;
  int id;
  boolean acertado;

  PontosInt InterReta(Reta r){
    return new PontosInt();
  };

  void Acertou(){
    this.acertado = true;
  }
}


// --------- Pixel ---------
class Pixel{
  RGB cor;
  PontosInt ptsInt;
  FigurasGeo PrimeiroObjeto;
  boolean visivel;

  Pixel(){
    this.cor = new RGB(30,144,255);
    this.visivel = true;
    this.ptsInt = null;
  }

  Pixel(PontosInt ptsInt){
    this.ptsInt = ptsInt;
    this.cor = ptsInt.pontos[0].cor;
    this.PrimeiroObjeto = ptsInt.PrimObj;
    this.visivel = true;
  }

  String printOrdem(){
    if (ptsInt != null){
      return ptsInt.ObjetosId();
    }
    else{
      return "Nao ha interseccoes";
    }
  }

  String printPontos(){
    if (ptsInt != null){
      return ptsInt.toString();
    }
    else{
      return "Nao ha interseccoes";
    }
  }
}


//----------- RGB ----------
class RGB{
  int R, G, B;
  RGB(int r, int g, int b){
    this.R = r;
    this.G = g;
    this.B = b;
  }
}


// ---------- Imagem ---------
class Imagem{
  Pixel pixels[][];
  int H, V;

  Imagem(int h, int v){
    this.pixels = new Pixel[h][v];
    for (int i = 0; i < h; i++){
      for (int j = 0; j < v; j++){
	pixels[i][j] = new Pixel();
      }
    }
    this.H = h;
    this.V = v;
  }

  void addPixel(Pixel p, int h, int v){
    if (h <= H && h >= 0 && v <= V && v >= 0){
      pixels[h][V-v-1] = p;
    }
    else{
      System.out.println("posicao na imagem("+h+", "+v+") invalida");
    }
  }

  void pixelNulo(int h, int v){
    if (h <= H && h >= 0 && v <= V && v >= 0){
      pixels[h][v] = new Pixel();
    }
    else{
      System.out.println("posicao na imagem("+h+", "+v+") invalida");
    }
  }

  Pixel getPixel(int h, int v){
    if (h <= H && h >= 0 && v <= V && v >= 0){
      return pixels[h][v];
    }
    else{
      System.out.println("posicao na imagem("+h+", "+v+") invalida");
    }
    return null;
  }

  void atingido(int h, int v){
    pixels[h][v].PrimeiroObjeto.Acertou();
    for (int i = 0; i < this.H; i++){
      for (int j = 0; j < this.V; j++){
        if (pixels[i][j].PrimeiroObjeto.acertado){
          pixels[i][j].visivel = true;
        }
      }
    }
  }
}

// ---------- Janela -------------
class Janela extends Frame{
  Janela(Espaco E, Muro M, Imagem I){

    JLabel l1 = new JLabel("Tiro:");
    JLabel l2 = new JLabel("Posicao do ponto no painel:");
    JLabel l3 = new JLabel("Ordem dos objetos atingidos:");
    JLabel l4 = new JLabel("Pontos atingidos:");
    JLabel posicao = new JLabel("");
    JTextArea ordem = new JTextArea("");
    JTextArea pontos = new JTextArea("");
    JLabel la = new JLabel("");

    l1.setBounds(500,40,100,30);
    l2.setBounds(500,110,250,30);
    posicao.setBounds(500,145,200,30);
    l3.setBounds(500,180,250,30);
    ordem.setBounds(500,215,250,30);
    l4.setBounds(500,250,250,30);
    pontos.setBounds(500,285,250,90);
    l1.setVisible(true);
    l2.setVisible(true);
    l3.setVisible(true);
    l4.setVisible(true);
    posicao.setVisible(true);
    ordem.setVisible(true);
    pontos.setVisible(true);

    JTextField TiroX = new JTextField();
    JTextField TiroY = new JTextField();

    TiroX.setBounds(500,80,50,20);
    TiroY.setBounds(550,80,50,20);

    JButton Atirar = new JButton("Atirar");
    JButton MostrarImagem = new JButton("Mostrar Imagem");

    Atirar.setBounds(620,70,100,30);
    MostrarImagem.setBounds(620,30,150,30);

    MostrarImagem.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e1) {
        for (int i = 0; i < I.H; i++){
          for (int j = 0; j < I.V; j++){
            I.atingido(i,j);
            Quadro aaa = new Quadro(I); 
          }
        }
      }
    });

    Atirar.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e2) {
        String s1 = TiroX.getText();
        String s2 = TiroY.getText();
        int h = Integer.parseInt(s1) - 1;
        int v = Integer.parseInt(s2) - 1;
        posicao.setText(M.buracos[h][v].toString());
        ordem.setText(I.pixels[h][v].printOrdem());
        pontos.setText(I.pixels[h][v].printPontos());
	I.atingido(h,v);
        Quadro aaa = new Quadro(I); 
      }
    });

    addWindowListener(new WindowAdapter(){
      public void windowClosing(WindowEvent f){
        dispose();  }    });

    add(l1); add(l2); add(l3); add(l4); add(posicao); add(ordem); add(pontos); add(TiroX); add(TiroY); add(Atirar); add(MostrarImagem); add(la);

    setVisible(true);
    setSize(900,500);
    setLocation(120,150);
  }
}

// ----------- Drawing Component --------------
class DrawingComponent extends Component{
  Imagem I;
  DrawingComponent(Imagem I){
    this.I = I;
  }
  public void paint(Graphics g){
    for (int i = 0; i < I.H; i++){
      for (int j = 0; j < I.V; j++){
        g.setColor(new Color(I.pixels[i][j].cor.R, I.pixels[i][j].cor.G, I.pixels[i][j].cor.B));
        g.drawRect(i,j,1,1);
      }
    }
  }
}
