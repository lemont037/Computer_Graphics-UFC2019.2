import java.awt.*;
import javax.swing.*;
import java.awt.event.*;
import java.lang.Math;
import java.text.DecimalFormat;
class Main {
  public static void main(String[] args) {
    Cenario C = new Cenario();
    Espaco E = new Espaco();
    E.AddObjeto(new Cilindro(new Ponto(0,-2,-10), new Vetor(0,1,0), 0.5, 2, new RGB(139,69,19), 1));
    E.AddObjeto(new Cone(new Ponto(0,0,-10), 3, 8, new Vetor(0,1,0), new RGB(0,255,127), 2));
    E.AddObjeto(new Cubo(new Ponto(0, 1,-20), 6, new RGB(178,34,34), 3));
    E.AddObjeto(new Cubo(new Ponto(0, 7,-20), 6, new RGB(178,34,34), 4));
    E.AddObjeto(new Cubo(new Ponto(0, 13,-20), 6, new RGB(178,34,34), 5));
    C.SetEspaco(E);
    Camera c = new Camera(new Ponto(0,0,0), new Ponto(0,0,-10),new Ponto(0,1,0));
    C.SetCamera(c);
    Janela janela = new Janela(C);
  }
}
class Camera{
  Ponto observador, LookAt, ViewUp;
  Vetor i, j, k;
  Muro FoV;

  Camera(Ponto Obs, Ponto LA, Ponto VU){
    this.observador = Obs;
    this.LookAt = LA;
    this.ViewUp = VU;
    Vetor K = observador.SubPonto(LA);
    this.k = K.ProdutoEscalar(1/Math.sqrt(K.ProdutoVetorial(K)));
    Vetor V = VU.SubPonto(Obs);
    Vetor I = new Vetor((V.Y*k.Z)-(V.Z*k.Y),(V.Z*k.X)-(V.X*k.Z),(V.X*k.Y)-(V.Y*k.X));
    this.i = I.ProdutoEscalar(1/Math.sqrt(I.ProdutoVetorial(I)));
    this.j = new Vetor((k.Y*i.Z)-(k.Z*i.Y), (k.Z*i.X)-(k.X*i.Z), (k.X*i.Y)-(k.Y*i.X));
    this.FoV = new Muro(4, k.SubPtVt(Obs), 500, 500);
  }

  Ponto MundoCamera(Ponto p){
    double Xc = (i.X*p.X) + (i.Y*p.Y) + (i.Z*p.Z) + (-observador.X*i.X);
    double Yc = (j.X*p.X) + (j.Y*p.Y) + (j.Z*p.Z) + (-observador.Y*j.Y);
    double Zc = (k.X*p.X) + (k.Y*p.Y) + (k.Z*p.Z) + (-observador.Z*k.Z);
    return new Ponto(Xc, Yc, Zc);
  }

  Vetor MundoCamera(Vetor v){
    double Xc = (i.X*v.X) + (i.Y*v.Y) + (i.Z*v.Z);
    double Yc = (j.X*v.X) + (j.Y*v.Y) + (j.Z*v.Z);
    double Zc = (k.X*v.X) + (k.Y*v.Y) + (k.Z*v.Z);
    return new Vetor(Xc, Yc, Zc);
  }

  Imagem CriarImagem(Espaco E){
    Espaco Ec = new Espaco();
    for (int k = 0; k < E.N; k++){
      Ec.AddObjeto(E.objetos[k].MundoCamera(this));
    }
    Imagem I = new Imagem(500,500);
    for (int i = 0; i < 500; i++){
      for (int j = 0; j < 500; j++){
        Ponto Aux = this.FoV.PnoMuro(i, j);
        Vetor N = Aux.SubPonto(this.observador);
        Vetor n = N.ProdutoEscalar(1/Math.sqrt(N.ProdutoVetorial(N)));
        Reta R = new Reta(this.observador, n);
        PontosInt Pfinal = new PontosInt();
        for (int l = 0; l < Ec.N; l++){
          PontosInt PInt = Ec.objetos[l].InterReta(R);
          Pfinal = Pfinal.Uniao(PInt);
        }
        if (Pfinal.N != 0){
          FigurasGeo prim = Ec.GetFigPorId(Pfinal.GetPonto(0).Id);
          Pfinal.SetPrimObj(prim);
          Pixel pix = new Pixel(Pfinal);
          I.addPixel(pix, i, j);
        }
        else{
          I.pixelNulo(i, j);
        }
      }
    }
    return I;
  }
}
class Cenario{
  Camera C;
  Espaco E;
  Imagem I;

  Cenario(Camera c, Espaco e, Imagem i){
    this.C = c;
    this.E = e;
    this.I = i;
  }

  Cenario(){
    this.E = new Espaco();
  }

  void SetCamera(Camera c){
    this.C = c;
  }

  void SetEspaco(Espaco e){
    this.E = e;
  }

  void SetImagem(Imagem i){
    this.I = i;
  }
}
class Cilindro extends FigurasGeo{
  
  Cilindro(Ponto B, Vetor u, double R, double H, RGB cor, int id){
    this.C = B;
    this.n = u;
    this.R = R;
    this.H = H;
    this.Cor = cor;
    this.Id = id;
    this.acertado = false;
  }

  double CalculoDeltaCil(Reta r){
    Vetor w = r.d.SubVetor(this.n.ProdutoEscalar(r.d.ProdutoVetorial(this.n))); // w = d - (d*u)*u
    Vetor v = r.P0.SubPonto(this.C).SubVetor(this.n.ProdutoEscalar(r.P0.SubPonto(this.C).ProdutoVetorial(this.n))); // v = (p0-B)-((p0-B)*u)*u
    double a = w.ProdutoVetorial(w); // a = (w*w)
    double b = v.ProdutoVetorial(w); // b = (v*w)
    double c = v.ProdutoVetorial(v) - (this.R * this.R); // c = (v*v) - R**2
    double Delta = (b*b) - (a*c);
    return Delta;
  }

  PontosInt InterReta(Reta r){ 
    Vetor w = r.d.SubVetor(this.n.ProdutoEscalar(r.d.ProdutoVetorial(this.n)));
    Vetor v = r.P0.SubPonto(this.C).SubVetor(this.n.ProdutoEscalar(r.P0.SubPonto(this.C).ProdutoVetorial(this.n)));
    double a = w.ProdutoVetorial(w);
    double b = v.ProdutoVetorial(w);
    double Delta = this.CalculoDeltaCil(r);
    if (Delta < 0){
      if (a == 0){
        Plano pl = new Plano(this.C, this.n);
        Ponto IntPlanoReta = pl.InterRetaPlano(r);
        double Dist = IntPlanoReta.Distancia(this.C);
        PontosInt PtsValidos = new PontosInt();
        if (Dist < this.R){
          Ponto topo = r.d.ProdutoEscalar(H).SomaPtVt(IntPlanoReta);
          topo.SetT(H + IntPlanoReta.T);
          this.SetAtributos(IntPlanoReta);
          this.SetAtributos(topo);
          PtsValidos.AddPonto(IntPlanoReta);
          PtsValidos.AddPonto(topo);
        }
        return PtsValidos;
      }
      return new PontosInt();
    }
    if (Delta == 0){
      double T = -b / a;
      Ponto p1 = r.PnaReta(T);
      double Verif = p1.SubPonto(this.C).ProdutoVetorial(this.n);
      if (Verif >= 0 && Verif <= H){
        p1.SetT(T);
        this.SetAtributos(p1);
        return new PontosInt(p1);
      }
    return new PontosInt();
    }
    if (Delta > 0){
      boolean p1V = false, p2V = false;
      double T = (-b + Math.sqrt(Delta)) / a;
      Ponto p1 = r.PnaReta(T);
      double Verif1 = p1.SubPonto(this.C).ProdutoVetorial(this.n);
      PontosInt PtsValidos = new PontosInt();
      if (Verif1 >= 0 && Verif1 <= H){
        p1.SetT(T);
        this.SetAtributos(p1);
        PtsValidos.AddPonto(p1);
        p1V = true;
      }
      T = (-b - Math.sqrt(Delta)) / a;
      Ponto p2 = r.PnaReta(T);
      double Verif2 = p2.SubPonto(this.C).ProdutoVetorial(this.n);
      if (Verif2 >= 0 && Verif2 <= H){
        p2.SetT(T);
        this.SetAtributos(p2);
        PtsValidos.AddPonto(p2);
        p2V = true;
      }
      if ((p1V && !(p2V)) || (!(p1V) && p2V)){ // Apenas um pt valido
        Plano pl = new Plano(this.C, this.n);
        Ponto IntPlReta = pl.InterRetaPlano(r);
        double Dist = IntPlReta.Distancia(C);
        if (Dist < this.R){
          this.SetAtributos(IntPlReta);
          PtsValidos.AddPonto(IntPlReta);
        }
        else{
          Ponto Int = r.d.ProdutoEscalar(H).SomaPtVt(IntPlReta);
          Int.SetT(H + IntPlReta.T);
          this.SetAtributos(Int);
          PtsValidos.AddPonto(Int);
        }
      }
      if (!(p1V) && !(p2V) && ((Verif1 < 0 && Verif2 > H)||(Verif1 > H && Verif2 < 0))){ // Dois invalidos e um acima e outro abaixo do cilindro
        Plano pl = new Plano(this.C, this.n);
        Ponto IntPlaReta = pl.InterRetaPlano(r);
        Ponto IntRetaTopo = r.d.ProdutoEscalar(H).SomaPtVt(IntPlaReta);
        IntRetaTopo.SetT(H + IntPlaReta.T);
        this.SetAtributos(IntPlaReta);
        this.SetAtributos(IntRetaTopo);
        PtsValidos.AddPonto(IntPlaReta);
        PtsValidos.AddPonto(IntRetaTopo);
      }
      return PtsValidos;
    }
    return new PontosInt();
  }

  Cilindro MundoCamera(Camera c){
    Ponto Cc = c.MundoCamera(this.C);
    Vetor uc = c.MundoCamera(this.n);
    return new Cilindro(Cc, uc, this.R, this.H, this.Cor, this.Id);
  }
}
class Cone extends FigurasGeo{ // cos**2(X) = H**2 / H**2 + R**2
  Ponto V; // V = vertice do cone

  Cone(Ponto C, Ponto V, double R, double H, Vetor n, RGB cor, int id){
    this.n = n;
    this.H = H;
    this.R = R;
    this.C = C;
    this.V = V;
    this.Cor = cor;
    this.Id = id;
    this.acertado = false;
  }

  Cone(Ponto C, double R, double H, Vetor n, RGB cor, int id){
    this.C = C;
    this.n = n;
    this.R = R;
    this.H = H;
    this.V = n.ProdutoEscalar(H).SomaPtVt(C);
    this.Cor = cor;
    this.Id = id;
    this.acertado = false;
  }

  Cone(double R, double H, Vetor n, Ponto V, RGB cor, int id){
    this.V = V;
    this.H = H;
    this.R = R;
    this.n = n;
    this.C = n.ProdutoEscalar(H).SubPtVt(V);
    this.Cor = cor;
    this.Id = id;
    this.acertado = false;
  }

  double CalculoDeltaCone(Reta r){
    Vetor v = V.SubPonto(r.P0);
    double a = (r.d.ProdutoVetorial(this.n))*(r.d.ProdutoVetorial(this.n)) - (r.d.ProdutoVetorial(r.d) * (H*H/((H*H)+(R*R)))); // (d*n)**2-(d*d)*cos**2(X)
    double b = v.ProdutoVetorial(r.d)*(H*H/((H*H)+(R*R)))-(v.ProdutoVetorial(n) * r.d.ProdutoVetorial(n)); // (v*d)*cos**2(X) - (v*n)*(d*n)
    double c = v.ProdutoVetorial(n)*v.ProdutoVetorial(n)-v.ProdutoVetorial(v)*(H*H/((H*H)+(R*R))); // (v*n)**2 - (v*v)*cos**2(X)
    double Delta = (b*b) - (a*c);
    return Delta;
  }

  PontosInt InterReta(Reta r){
    Vetor v = V.SubPonto(r.P0);
    double a = (r.d.ProdutoVetorial(this.n))*(r.d.ProdutoVetorial(this.n)) - (r.d.ProdutoVetorial(r.d) * (H*H/((H*H)+(R*R))));
    double b = v.ProdutoVetorial(r.d)*(H*H/((H*H)+(R*R)))-(v.ProdutoVetorial(n) * r.d.ProdutoVetorial(n));
    double c = v.ProdutoVetorial(n)*v.ProdutoVetorial(n)-v.ProdutoVetorial(v)*(H*H/((H*H)+(R*R)));
    double Delta = this.CalculoDeltaCone(r);
    if (a == 0){
      double T = -(c / (2 * b));
      Ponto p1 = r.PnaReta(T);
      if (V.SubPonto(p1).ProdutoVetorial(n) >= 0 && V.SubPonto(p1).ProdutoVetorial(n) <= H){
        p1.SetT(T);
        this.SetAtributos(p1);
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
      if (V.SubPonto(p1).ProdutoVetorial(n) >= 0 && V.SubPonto(p1).ProdutoVetorial(n) <= H){
        p1.SetT(T);
        this.SetAtributos(p1);
        return new PontosInt(p1);
      }
      return null;
    }
    if (Delta > 0){
      PontosInt IntRetaCon = new PontosInt();
      boolean PV1 = false, PV2 = false;
      double T1 = (-b + Math.sqrt(Delta)) / a;
      Ponto p1 = r.PnaReta(T1);
      if (V.SubPonto(p1).ProdutoVetorial(n) >= 0 && V.SubPonto(p1).ProdutoVetorial(n) <= H){
        PV1 = true;
        p1.SetT(T1);
        this.SetAtributos(p1);
        IntRetaCon.AddPonto(p1);
      }
      double T2 = (-b - Math.sqrt(Delta)) / a;
      Ponto p2 = r.PnaReta(T2);
      if (V.SubPonto(p2).ProdutoVetorial(n) >= 0 && V.SubPonto(p2).ProdutoVetorial(n) <= H){
        PV2 = true;
        p2.SetT(T2);
        this.SetAtributos(p2);
        IntRetaCon.AddPonto(p2);
      }
      Plano pl = new Plano(this.C, this.n);
      Ponto IntRetaPl = pl.InterRetaPlano(r);
      if(IntRetaPl.Distancia(this.C) < this.R){
        this.SetAtributos(IntRetaPl);
        IntRetaCon.AddPonto(IntRetaPl);
      }
      return IntRetaCon;
    }
    return new PontosInt();
  }

  Cone MundoCamera(Camera c){
    Ponto Cc = c.MundoCamera(this.C);
    Vetor nc = c.MundoCamera(this.n);
    return new Cone(Cc, this.R, this.H, nc, this.Cor, this.Id);
  }
}
class Cubo extends FigurasGeo{
  double A; // aresta
  Ponto p1, p2, p3, p4, p5, p6, p7, p8;

  Cubo(Ponto p, double A, RGB cor, int id){
    this.C = p;
    this.A = A;
    this.Cor = cor;
    this.Id = id;
    this.acertado = false;
    this.n = null;
    this.R = A;
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
    Ponto P1 = new Ponto(C.X + A/2, C.Y, C.Z);
    Ponto P2 = new Ponto(C.X - A/2, C.Y, C.Z);
    Ponto P3 = new Ponto(C.X, C.Y + A/2, C.Z);
    Ponto P4 = new Ponto(C.X, C.Y - A/2, C.Z);
    Ponto P5 = new Ponto(C.X, C.Y, C.Z + A/2);
    Ponto P6 = new Ponto(C.X, C.Y, C.Z - A/2);
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
      this.SetAtributos(Int1);
      ptsInt.AddPonto(Int1);
    }
    if (Math.abs(Int2.X - P2.X) <= A/2 && Math.abs(Int2.Y - P2.Y) <= A/2 && Math.abs(Int2.Z - P2.Z) <= A/2){
      this.SetAtributos(Int2);
      ptsInt.AddPonto(Int2);
    }
    if (Math.abs(Int3.X - P3.X) <= A/2 && Math.abs(Int3.Y - P3.Y) <= A/2 && Math.abs(Int3.Z - P3.Z) <= A/2){
      this.SetAtributos(Int3);
      ptsInt.AddPonto(Int3);
    }
    if (Math.abs(Int4.X - P4.X) <= A/2 && Math.abs(Int4.Y - P4.Y) <= A/2 && Math.abs(Int4.Z - P4.Z) <= A/2){
      this.SetAtributos(Int4);
      ptsInt.AddPonto(Int4);
    }
    if (Math.abs(Int5.X - P5.X) <= A/2 && Math.abs(Int5.Y - P5.Y) <= A/2 && Math.abs(Int5.Z - P5.Z) <= A/2){
      this.SetAtributos(Int5);
      ptsInt.AddPonto(Int5);
    }
    if (Math.abs(Int6.X - P6.X) <= A/2 && Math.abs(Int6.Y - P6.Y) <= A/2 && Math.abs(Int6.Z - P6.Z) <= A/2){
      this.SetAtributos(Int6);
      ptsInt.AddPonto(Int6);
    }
    return ptsInt;
  }

  Cubo MundoCamera(Camera c){
    Ponto Pc = c.MundoCamera(this.C);
    return new Cubo(Pc, this.A, this.Cor, this.Id);
  }
}
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
class Esfera extends FigurasGeo{

  Esfera(Ponto C, double R, RGB cor, int id){
    this.C = C;
    this.R = R;
    this.Cor = cor;
    this.Id = id;
    this.acertado = false;
  }

  double CalculoDeltaEsf(Reta R){
    double a = R.d.ProdutoVetorial(R.d); // a = (d*d)
    double b = (R.P0.SubPonto(this.C)).ProdutoVetorial(R.d); // b = ((p0-C)*d)
    double c = (R.P0.SubPonto(this.C)).ProdutoVetorial(R.P0.SubPonto(this.C)) - (this.R * this.R); // c = ((p0-C)*(p0-C)-(R**2))
    double Delta = (b*b) - (a*c);
    return Delta;
  }

  PontosInt InterReta(Reta r){ // (P-C)*(P-C) = R**2
    double a = r.d.ProdutoVetorial(r.d);
    double b = (r.P0.SubPonto(this.C)).ProdutoVetorial(r.d);
    double Delta = this.CalculoDeltaEsf(r);
    if (Delta < 0) {
      new PontosInt();
    }
    if (Delta == 0){
      double T = -b/a;
      Ponto p1 = r.PnaReta(T);
      p1.SetT(T);
      this.SetAtributos(p1);
      return new PontosInt(p1);
    }
    if (Delta > 0){
      double T = (-b + Math.sqrt(Delta)) / a;
      Ponto p1 = r.PnaReta(T);
      p1.SetT(T);
      this.SetAtributos(p1);
      T = (-b - Math.sqrt(Delta)) / a;
      Ponto p2 = r.PnaReta(T);
      p2.SetT(T);
      this.SetAtributos(p2);
      return new PontosInt(p1, p2);
    }
    return new PontosInt();
  }

  Esfera MundoCamera(Camera c){
    Ponto Cc = c.MundoCamera(this.C);
    return new Esfera(Cc, this.R, this.Cor, this.Id);
  }
}
class Espaco{
  int N; // N = numero de objetos no Espaco
  FigurasGeo objetos[];

  Espaco(){
    this.N = 0;
    this.objetos = new FigurasGeo[10];
  }

  void AddObjeto(FigurasGeo obj){
    objetos[this.N] = obj;
    this.N += 1;
  }

  FigurasGeo GetFigPorId(int id){
    for (int i = 0; i < N; i++){
      if (objetos[i].Id == id){
        return objetos[i];
      }
    }
    return null;
  }
}
abstract class FigurasGeo{
  RGB Cor;
  int Id;
  boolean acertado;
  Ponto C;
  double R, H;
  Vetor n;

  PontosInt InterReta(Reta r){
    return new PontosInt();
  }

  void Acertou(){
    this.acertado = true;
  }

  void SetAtributos(Ponto p){
    p.SetCor(this.Cor);
    p.SetId(this.Id);
  }

  FigurasGeo MundoCamera(Camera c){
    return null;
  }
}
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
class Janela extends Frame{
  Janela(Cenario C){

    JLabel Observador = new JLabel("Observador:");
    JLabel LookAt = new JLabel("Look At:");
    JLabel ViewUp = new JLabel("View Up:");
    JLabel la = new JLabel("");

    JTextField ObsX = new JTextField("");
    JTextField ObsY = new JTextField("");
    JTextField ObsZ = new JTextField("");
    JTextField LAX = new JTextField("");
    JTextField LAY = new JTextField("");
    JTextField LAZ = new JTextField("");
    JTextField VUX = new JTextField("");
    JTextField VUY = new JTextField("");
    JTextField VUZ = new JTextField("");

    JButton MostrarImagem = new JButton("Mostrar Imagem");
    JButton CriaCone = new JButton("Criar Cone");
    JButton CriaCilindro = new JButton("Criar Cilindro");
    JButton CriaEsfera = new JButton("Criar Esfera");
    JButton CriaCubo = new JButton("Criar Cubo");
    JButton PosicionarObs = new JButton("Posicionar Observador");

    Observador.setBounds(30, 30, 100, 30);
    LookAt.setBounds(30, 90, 100, 30);
    ViewUp.setBounds(30, 150, 100, 30);
    ObsX.setBounds(30, 60, 30, 20);
    ObsY.setBounds(60, 60, 30, 20);
    ObsZ.setBounds(90, 60, 30, 20);
    LAX.setBounds(30, 120, 30, 20);
    LAY.setBounds(60, 120, 30, 20);
    LAZ.setBounds(90, 120, 30, 20);
    VUX.setBounds(30, 180, 30, 20);
    VUY.setBounds(60, 180, 30, 20);
    VUZ.setBounds(90, 180, 30, 20);
    MostrarImagem.setBounds(30, 250, 150, 30);
    CriaCone.setBounds(30, 290, 150, 30);
    CriaCilindro.setBounds(190, 290, 150, 30);
    CriaCubo.setBounds(30, 325, 150, 30);
    CriaEsfera.setBounds(190, 325, 150, 30);
    PosicionarObs.setBounds(30, 210, 200, 30);

    add(Observador);add(LookAt);add(ViewUp);add(ObsX);add(ObsY);add(ObsZ);add(LAX);add(LAY);add(LAZ);add(VUX);add(VUY);add(VUZ);add(MostrarImagem);add(CriaCone);add(CriaCilindro);add(CriaEsfera);add(CriaCubo);add(PosicionarObs);add(la);

    MostrarImagem.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e1){
        try{
          C.SetImagem(C.C.CriarImagem(C.E));
          Quadro Q = new Quadro(C);
        }
        catch(Exception f4){}
      }
    });

    CriaCone.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e2){
        JanelaCCone aux = new JanelaCCone(C);
      }
    });

    CriaCilindro.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e2){
        JanelaCCilindro aux = new JanelaCCilindro(C);
      }
    });

    CriaCubo.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e2){
        JanelaCCubo aux = new JanelaCCubo(C);
      }
    });

    CriaEsfera.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e2){
        JanelaCEsfera aux = new JanelaCEsfera(C);
      }
    });

    PosicionarObs.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e2){
        String obsx = ObsX.getText();
        String obsy = ObsY.getText();
        String obsz = ObsZ.getText();
        String lax = LAX.getText();
        String lay = LAY.getText();
        String laz = LAZ.getText();
        String vux = VUX.getText();
        String vuy = VUY.getText();
        String vuz = VUZ.getText();
        try{
          double obsX = Double.parseDouble(obsx);
          double obsY = Double.parseDouble(obsy);
          double obsZ = Double.parseDouble(obsz);
          double laX = Double.parseDouble(lax);
          double laY = Double.parseDouble(lay);
          double laZ = Double.parseDouble(laz);
          double vuX = Double.parseDouble(vux);
          double vuY = Double.parseDouble(vuy);
          double vuZ = Double.parseDouble(vuz);
          Ponto observador = new Ponto(obsX, obsY, obsZ);
          Ponto lookat = new Ponto(laX, laY, laZ);
          Ponto viewup = new Ponto(vuX, vuY, vuZ);
          Camera c = new Camera(observador, lookat, viewup);
          C.SetCamera(c);
        }
        catch(Exception f){ //talvez colocar mensagem de erro
        }
      }
    });

    addWindowListener(new WindowAdapter(){
      public void windowClosing(WindowEvent f){
        dispose();  }    });

    setVisible(true);
    setSize(400,400);
    setLocation(120,150);
  }
}
class JanelaCCilindro extends Frame{
  JanelaCCilindro(Cenario c){
    JLabel ponto = new JLabel("Ponto do Centro:");
    JLabel vetor = new JLabel("Vetor de Direcao:");
    JLabel raio = new JLabel("Raio:");
    JLabel altura = new JLabel("Altura:");
    JLabel id = new JLabel("Id:");
    JLabel reflexividade = new JLabel("Reflexividade:");
    JTextField pontox = new JTextField("");
    JTextField pontoy = new JTextField("");
    JTextField pontoz = new JTextField("");
    JTextField vetorx = new JTextField("");
    JTextField vetory = new JTextField("");
    JTextField vetorz = new JTextField("");
    JTextField r = new JTextField("");
    JTextField h = new JTextField("");
    JTextField Id = new JTextField("");
    JTextField refleX = new JTextField("");
    JTextField refleY = new JTextField("");
    JTextField refleZ = new JTextField("");
    JButton CriaCilindro = new JButton("Criar Cilindro");

    /*ponto.setBounds(x, y, width, height);
    vetor.setBounds(x, y, width, height);
    raio.setBounds(x, y, width, height);
    altura.setBounds(x, y, width, height);
    id.setBounds(x, y, width, height);
    reflexividade.setBounds(x, y, width, height);
    pontox.setBounds(x, y, width, height);
    pontoy.setBounds(x, y, width, height);
    pontoz.setBounds(x, y, width, height);
    vetorx.setBounds(x, y, width, height);
    vetory.setBounds(x, y, width, height);
    vetorz.setBounds(x, y, width, height);
    r.setBounds(x, y, width, height);
    h.setBounds(x, y, width, height);
    Id.setBounds(x, y, width, height);
    refleX.setBounds(x, y, width, height);
    refleY.setBounds(x, y, width, height);
    refleZ.setBounds(x, y, width, height);
    CriaCilindro.setBounds(x, y, width, height);*/

    CriaCilindro.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
        String px = pontox.getText();
        String py = pontoy.getText();
        String pz = pontoz.getText();
        String vx = vetorx.getText();
        String vy = vetory.getText();
        String vz = vetorz.getText();
        String Raio = r.getText();
        String Altura = h.getText();
        String ID = Id.getText();
        String reflex = refleX.getText();
        String refley = refleY.getText();
        String reflez = refleZ.getText();
        try{
          double Px = Double.parseDouble(px);
          double Py = Double.parseDouble(py);
          double Pz = Double.parseDouble(pz);
          double Vx = Double.parseDouble(vx);
          double Vy = Double.parseDouble(vy);
          double Vz = Double.parseDouble(vz);
          double R = Double.parseDouble(Raio);
          double H = Double.parseDouble(Altura);
          int I = Integer.parseInt(ID);
          int RefX = Integer.parseInt(reflex);
          int RefY = Integer.parseInt(refley);
          int RefZ = Integer.parseInt(reflez);
          RGB cor = new RGB(RefX, RefY, RefZ);
          c.E.AddObjeto(new Cilindro(new Ponto(Px, Py, Pz), new Vetor(Vx, Vy, Vz), R, H, cor, I));
          dispose();
        }
        catch(Exception f){
        }
      }
    });
  
    addWindowListener(new WindowAdapter(){
      public void windowClosing(WindowEvent f){
        dispose();  }    });
  
    setVisible(true);
    setSize(500,500);
    setLocation(120,150);
  }
}
class JanelaCCone extends Frame{
  JanelaCCone(Cenario c){
    JLabel ponto = new JLabel("Ponto do Centro:");
    JLabel vetor = new JLabel("Vetor de Direcao:");
    JLabel raio = new JLabel("Raio:");
    JLabel altura = new JLabel("Altura:");
    JLabel id = new JLabel("Id:");
    JLabel reflexividade = new JLabel("Reflexividade:");
    JTextField pontox = new JTextField("");
    JTextField pontoy = new JTextField("");
    JTextField pontoz = new JTextField("");
    JTextField vetorx = new JTextField("");
    JTextField vetory = new JTextField("");
    JTextField vetorz = new JTextField("");
    JTextField r = new JTextField("");
    JTextField h = new JTextField("");
    JTextField Id = new JTextField("");
    JTextField refleX = new JTextField("");
    JTextField refleY = new JTextField("");
    JTextField refleZ = new JTextField("");
    JButton CriaCone = new JButton("Criar Cone");

    CriaCone.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
        String px = pontox.getText();
        String py = pontoy.getText();
        String pz = pontoz.getText();
        String vx = vetorx.getText();
        String vy = vetory.getText();
        String vz = vetorz.getText();
        String Raio = r.getText();
        String Altura = h.getText();
        String ID = Id.getText();
        String reflex = refleX.getText();
        String refley = refleY.getText();
        String reflez = refleZ.getText();
        try{
          double Px = Double.parseDouble(px);
          double Py = Double.parseDouble(py);
          double Pz = Double.parseDouble(pz);
          double Vx = Double.parseDouble(vx);
          double Vy = Double.parseDouble(vy);
          double Vz = Double.parseDouble(vz);
          double R = Double.parseDouble(Raio);
          double H = Double.parseDouble(Altura);
          int I = Integer.parseInt(ID);
          int RefX = Integer.parseInt(reflex);
          int RefY = Integer.parseInt(refley);
          int RefZ = Integer.parseInt(reflez);
          RGB cor = new RGB(RefX, RefY, RefZ);
          c.E.AddObjeto(new Cone(new Ponto(Px, Py, Pz), R, H, new Vetor(Vx, Vy, Vz), cor, I));
          dispose();
        }
        catch(Exception f){
        }
      }
    });

    addWindowListener(new WindowAdapter(){
      public void windowClosing(WindowEvent f){
        dispose();  }    });
    
    setVisible(true);
    setSize(500,500);
    setLocation(120,150);
  }
}
class JanelaCCubo extends Frame{
  JanelaCCubo(Cenario c){
    JLabel ponto = new JLabel("Ponto do Centro:");
    JLabel lado = new JLabel("Lado:");
    JLabel id = new JLabel("Id:");
    JLabel reflexividade = new JLabel("Reflexividade:");
    JTextField pontox = new JTextField("");
    JTextField pontoy = new JTextField("");
    JTextField pontoz = new JTextField("");
    JTextField l = new JTextField("");
    JTextField Id = new JTextField("");
    JTextField refleX = new JTextField("");
    JTextField refleY = new JTextField("");
    JTextField refleZ = new JTextField("");
    JButton CriaCubo = new JButton("Criar Cubo");

    CriaCubo.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
        String px = pontox.getText();
        String py = pontoy.getText();
        String pz = pontoz.getText();
        String Lado = l.getText();
        String ID = Id.getText();
        String reflex = refleX.getText();
        String refley = refleY.getText();
        String reflez = refleZ.getText();
        try{
          double Px = Double.parseDouble(px);
          double Py = Double.parseDouble(py);
          double Pz = Double.parseDouble(pz);
          double L = Double.parseDouble(Lado);
          int I = Integer.parseInt(ID);
          int RefX = Integer.parseInt(reflex);
          int RefY = Integer.parseInt(refley);
          int RefZ = Integer.parseInt(reflez);
          RGB cor = new RGB(RefX, RefY, RefZ);
          c.E.AddObjeto(new Cubo(new Ponto(Px, Py, Pz), L, cor, I));
        }
        catch(Exception f){
        }
      }
    });

    addWindowListener(new WindowAdapter(){
      public void windowClosing(WindowEvent f){
        dispose();  }    });
    
    setVisible(true);
    setSize(500,500);
    setLocation(120,150);
  }
}
class JanelaCEsfera extends Frame{
  JanelaCEsfera(Cenario c){
    JLabel ponto = new JLabel("Ponto do Centro:");
    JLabel raio = new JLabel("Raio:");
    JLabel id = new JLabel("Id:");
    JLabel reflexividade = new JLabel("Reflexividade:");
    JTextField pontox = new JTextField("");
    JTextField pontoy = new JTextField("");
    JTextField pontoz = new JTextField("");
    JTextField r = new JTextField("");
    JTextField Id = new JTextField("");
    JTextField refleX = new JTextField("");
    JTextField refleY = new JTextField("");
    JTextField refleZ = new JTextField("");
    JButton CriaEsfera = new JButton("Criar Esfera");

    CriaEsfera.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e){
        String px = pontox.getText();
        String py = pontoy.getText();
        String pz = pontoz.getText();
        String Raio = r.getText();
        String ID = Id.getText();
        String reflex = refleX.getText();
        String refley = refleY.getText();
        String reflez = refleZ.getText();
        try{
          double Px = Double.parseDouble(px);
          double Py = Double.parseDouble(py);
          double Pz = Double.parseDouble(pz);
          double R = Double.parseDouble(Raio);
          int I = Integer.parseInt(ID);
          int RefX = Integer.parseInt(reflex);
          int RefY = Integer.parseInt(refley);
          int RefZ = Integer.parseInt(reflez);
          RGB cor = new RGB(RefX, RefY, RefZ);
          c.E.AddObjeto(new Esfera(new Ponto(Px, Py, Pz), R, cor, I));
          dispose();
        }
        catch(Exception f){
        }
      }
    });

    addWindowListener(new WindowAdapter(){
      public void windowClosing(WindowEvent f){
        dispose();  }    });

    setVisible(true);
    setSize(500,500);
    setLocation(120,150);
  }
}
class JanelaObj extends Frame{
  JanelaObj(FigurasGeo O){
    JLabel centro = new JLabel("Centro");
    JTextField Centro = new JTextField("");
    JLabel raio = new JLabel("Raio");
    JTextField Raio = new JTextField("");
    JLabel altura = new JLabel("Altura");
    JTextField Altura = new JTextField("");
    JLabel VDir = new JLabel("Vetor de Direcao");
    JTextField VetorDir = new JTextField("");
    JLabel Reflexividade = new JLabel("Reflexividade");
    JTextField reflexividade = new JTextField("");
    JLabel id = new JLabel("Id");
    JTextField Id = new JTextField("");
    JLabel lado = new JLabel("Lado");
    JLabel la = new JLabel("");
    centro.setVisible(true);
    Centro.setVisible(true);
    raio.setVisible(false);
    Raio.setVisible(false);
    altura.setVisible(false);
    Altura.setVisible(false);
    VDir.setVisible(false);
    VetorDir.setVisible(false);
    Reflexividade.setVisible(true);
    reflexividade.setVisible(true);
    id.setVisible(true);
    Id.setVisible(true);
    lado.setVisible(false);
    reflexividade.setText(O.Cor.toString());
    Id.setText(Integer.toString(O.Id));
    Centro.setText(O.C.toString());
    centro.setBounds(30, 30, 100, 30);
    Centro.setBounds(30, 60, 150, 30);
    raio.setBounds(200, 100, 100, 30);
    Raio.setBounds(200, 130, 50, 30);
    altura.setBounds(200, 30, 100, 30);
    Altura.setBounds(200, 30, 50, 30);
    VDir.setBounds(30, 100, 150, 30);
    VetorDir.setBounds(30, 130, 150, 30);
    Reflexividade.setBounds(30, 170, 100, 30);
    reflexividade.setBounds(30, 210, 150, 30);
    id.setBounds(200, 170, 50, 30);
    Id.setBounds(200, 210, 50, 30);
    lado.setBounds(200, 100, 100, 30);
    add(Centro);add(centro);add(raio);add(Raio);add(altura);add(Altura);add(VDir);add(VetorDir);add(Reflexividade);add(reflexividade);add(id);add(Id);add(lado);add(la);

    if (O instanceof Cone){
      raio.setVisible(true);
      Raio.setVisible(true);
      Raio.setText(Double.toString(O.R));
      altura.setVisible(true);
      Altura.setVisible(true);
      Altura.setText(Double.toString(O.H));
      VDir.setVisible(true);
      VetorDir.setVisible(true);
      VetorDir.setText(O.n.toString());
    }
    if (O instanceof Cilindro){
      raio.setVisible(true);
      Raio.setVisible(true);
      Raio.setText(Double.toString(O.R));
      altura.setVisible(true);
      Altura.setVisible(true);
      Altura.setText(Double.toString(O.H));
      VDir.setVisible(true);
      VetorDir.setVisible(true);
      VetorDir.setText(O.n.toString());
    }
    if (O instanceof Esfera){
      raio.setVisible(true);
      Raio.setVisible(true);
      Raio.setText(Double.toString(O.R));
    }
    if (O instanceof Cubo){
      lado.setVisible(true);
      Raio.setVisible(true);
      Raio.setText(Double.toString(O.R));
    }
    setVisible(true);
    setSize(400,400);
    setLocation(120,150);
  }
}
class Muro{
  double L; // L = lado
  Ponto P;
  Ponto buracos[][];
  int H, V;

  Muro(double L, Ponto C, int h, int v){
    this.P = new Ponto((C.X - L/2), (C.Y + L/2), C.Z);
    this.L = L;
    this.buracos = new Ponto[h][v];
    this.H = h;
    this.V = v;
  }

  Ponto PnoMuro(int h, int v){ // Calculo das coords X e Y do ponto no muro, L/2H+(h-1)*L/H
    Ponto p = new Ponto(this.P.X + (this.L/(2*this.H)) + (h*this.L/this.H), this.P.Y - (this.L/(2*this.V)) - (v*this.L/this.V), this.P.Z);
    this.buracos[h][v] = p;
    return p;
  }
}
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
    this.cor = ptsInt.pontos[0].Cor;
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
class Plano{
  Ponto PPl; // Ponto qualquer do plano
  Vetor n; // Vetor unitario perpendicular ao plano

  Plano(Ponto P, Vetor n){
    this.PPl = P;
    this.n = n;
  }

  Ponto InterRetaPlano(Reta R){ // (P - Ppl) * n = 0
    if (this.n.ProdutoVetorial(R.d) == 0){
      return null;
    }
    Vetor PPlP0 = this.PPl.SubPonto(R.P0);
    double T =  PPlP0.ProdutoVetorial(this.n) / this.n.ProdutoVetorial(R.d);
    Ponto PInt = R.PnaReta(T);
    PInt.SetT(T);
    return PInt;
  }
}
class Ponto{
  double X, Y, Z, T; // Coordenadas X,Y,Z e caso seja um pt de intersecao escalar T
  int Id; // id do objeto caso seja um pt de intersecao
  RGB Cor; // cor do objeto caso seja um pt de intersecao

  Ponto(double x, double y, double z){
    this.X = x;
    this.Y = y;
    this.Z = z;
  }

  void SetPonto(double x, double y, double z){
    this.X = x;
    this.Y = y;
    this.Z = z;
  }

  void SetT(double t){
    this.T = t;
  }

  void SetId(int id){
    this.Id = id;
  }

  void SetCor(RGB cor){
    this.Cor = cor;
  }

  Ponto SomaPonto(Ponto p){
    double XSoma = this.X + p.X;
    double YSoma = this.Y + p.Y;
    double ZSoma = this.Z + p.Z;
    return new Ponto(XSoma, YSoma, ZSoma);
  }

  Vetor SubPonto(Ponto p){
    double XSub = this.X - p.X;
    double YSub = this.Y - p.Y;
    double ZSub = this.Z - p.Z;
    return new Vetor(XSub, YSub, ZSub);
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

  String ID(){
    return Integer.toString(this.Id);
  }
}
class PontosInt{ // Classe feita para guardar os pontos de intersecao
  Ponto pontos[];
  int N; // numero de pts
  FigurasGeo PrimObj; // primeiro objeto que os pontos interceptam

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

  Ponto GetPonto(int posicao){
    if (this.N > 0){
      return this.pontos[posicao];
    }
    return null;
  }

  void SetPrimObj(FigurasGeo p){
    this.PrimObj = p;
  }

  void AddPonto(Ponto p){
    if (this.N == 0){
      this.pontos[0] = p;
      this.N = 1;
    }
    else{
      int i = 0;
      while (i < this.N){
        if (p.T >= this.pontos[i].T){
          int j = this.N;
          while (j > i){
            this.pontos[j] = this.pontos[j-1];
            j = j - 1;
          }
          this.pontos[j] = p;
          i = this.N;
        }
        i += 1;
      }
      if (i == this.N){
        this.pontos[i] = p;
      }
      this.N += 1;
    }
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
      uniao.AddPonto(p.GetPonto(i));
    }
    for (int j = 0; j < this.N; j++){
      uniao.AddPonto(this.GetPonto(j));
    }
    return uniao;
  }

  String ObjetosId(){
    if (this.N == 0){
      return "Nao ha intersecao";
    }
    String Ids = pontos[0].ID();
    for (int i = 1; i < this.N; i++){
      String aux1 = pontos[i-1].ID();
      String aux2 = pontos[i].ID();
      if (!(aux1.equals(aux2))){
        Ids = Ids.concat(" " + aux2);
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
class Quadro extends Frame{
  Quadro(Cenario C){
    add(new DrawingComponent(C.I));
    addMouseListener(new MouseAdapter(){ // terminar funcao do click
      public void mouseCliked(MouseEvent e){
        int X = e.getX();
        int Y = e.getY();
	System.out.println("" + X + ", " + Y);
        try{
          JanelaObj aux = new JanelaObj(C.I.pixels[X][Y].PrimeiroObjeto);
        }
        catch(Exception e1){}
      }
    });
    setVisible(true);
    setSize(C.I.H,C.I.V);
    setLocation(520,150);
    addWindowListener(new WindowAdapter(){
      public void windowClosing(WindowEvent f){
        dispose();  }    });
  }
}
class Reta{
  Ponto P0; // Ponto inicial da reta
  Vetor d; // Vetor unitario de sentido e direcao da reta

  Reta(Ponto p0, Vetor d){
    this.P0 = p0;
    this.d = d;
  }

  Ponto PnaReta(double t){ // P(t) = p0 + t*d
    Vetor V = this.d.ProdutoEscalar(t);
    Ponto P =  V.SomaPtVt(this.P0);
    return P;
  }
}
class RGB{
  int R, G, B;
  RGB(int r, int g, int b){
    this.R = r;
    this.G = g;
    this.B = b;
  }

  public String toString(){
    return "("+this.R+", "+this.G+", "+this.B+")";
  }
}
class Vetor{
  double X, Y, Z;

  Vetor(double x, double y, double z){
    this.X = x;
    this.Y = y;
    this.Z = z;
  }

  Vetor ProdutoEscalar(double esc){
    double XEsc = this.X * esc;
    double YEsc = this.Y * esc;
    double ZEsc = this.Z * esc;
    return new Vetor(XEsc, YEsc, ZEsc);
  }

  Vetor ProdutoEscalar(Vetor v){
    double XEsc = this.X * v.X;
    double YEsc = this.Y * v.Y;
    double ZEsc = this.Z * v.Z;
    return new Vetor(XEsc, YEsc, ZEsc);
  }

  double ProdutoVetorial(Vetor v){
    double Prod = (this.X*v.X) + (this.Y*v.Y) + (this.Z*v.Z);
    return Prod;
  }

  double ProdutoVetorial(Ponto p){
    double Prod = (this.X*p.X) + (this.Y*p.Y) + (this.Z*p.Z);
    return Prod;
  }

  Ponto SomaPtVt(Ponto p){
    double XSoma = this.X + p.X;
    double YSoma = this.Y + p.Y;
    double ZSoma = this.Z + p.Z;
    return new Ponto(XSoma, YSoma, ZSoma);
  }

  Ponto SubPtVt(Ponto p){
    double XSub = p.X - this.X;
    double YSub = p.Y - this.Y;
    double ZSub = p.Z - this.Z;
    return new Ponto(XSub, YSub, ZSub);
  }

  Vetor SubVetor(Vetor v){
    double XSub = this.X - v.X;
    double YSub = this.Y - v.Y;
    double ZSub = this.Z - v.Z;
    return new Vetor(XSub, YSub, ZSub);
  }

  public String toString(){
    return "("+this.X+", "+this.Y+", "+this.Z+")";
  }
}
