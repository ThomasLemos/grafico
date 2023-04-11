#include "stdio.h"
#include "raylib.h"
#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include <stdlib.h>
#include <math.h>
//aperte F5 para compilar
// Nome do objeto a ser carregado
char filename[100] = "objetos/maca2.byu";
// Aperte r para aplicar a visão da camera

const int screenWidth = 800;
const int screenHeight = 450;

float obj[3][3][10000];
float camera[3][14];
float Phong[3][3][10000];
int zbuffer[800][450];

float AreaTriangle(float P_1[2], float P_2[2], float P_3[2])
{
    return (((P_1[0] - P_2[0]) * (P_1[1] - P_3[1])) - ((P_1[0] - P_3[0]) * (P_1[1] - P_2[1]))) / 2;
}

void computePoint(float P1[3], float P2[3], float P3[3], float a, float b, float c, float pontoOriginal[2][3], int d)
{
    pontoOriginal[d][0] = (P1[0] * a + P2[0] * b + P3[0] * c);
    pontoOriginal[d][1] = (P1[1] * a + P2[1] * b + P3[1] * c);
    pontoOriginal[d][2] = (P1[2] * a + P2[2] * b + P3[2] * c);
}

void pOriginal(int V[3], float P_[2], float P_1[2], float P_2[2], float P_3[2], float pontoOriginal[2][3])
{
    float AreaT = AreaTriangle(P_1, P_2, P_3);
    float AreaA = AreaTriangle(P_, P_2, P_3);
    float AreaB = AreaTriangle(P_1, P_, P_3);
    float AreaC = AreaTriangle(P_1, P_2, P_);

    float a = AreaA / AreaT;
    float b = AreaB / AreaT;
    float c = AreaC / AreaT;

    float P1[3] = {Phong[0][0][V[0]], Phong[0][1][V[0]], Phong[0][2][V[0]]};
    float P2[3] = {Phong[0][0][V[1]], Phong[0][1][V[1]], Phong[0][2][V[1]]};
    float P3[3] = {Phong[0][0][V[2]], Phong[0][1][V[2]], Phong[0][2][V[2]]};

    float normalP1[3] = {Phong[2][0][V[0]], Phong[2][1][V[0]], Phong[2][2][V[0]]};
    float normalP2[3] = {Phong[2][0][V[1]], Phong[2][1][V[1]], Phong[2][2][V[1]]};
    float normalP3[3] = {Phong[2][0][V[2]], Phong[2][1][V[2]], Phong[2][2][V[2]]};

    computePoint(P1, P2, P3, a, b, c, pontoOriginal, 0);
    computePoint(normalP1, normalP2, normalP3, a, b, c, pontoOriginal, 1);
}

void normalize(float V[3], float C)
{
    V[0] = (1 / C) * V[0];
    V[1] = (1 / C) * V[1];
    V[2] = (1 / C) * V[2];
}

void normalize2(float V1[3], float F, float V2[3])
{
    V1[0] = (1 / F) * V2[0];
    V1[1] = (1 / F) * V2[1];
    V1[2] = (1 / F) * V2[2];
}

void normalize3(float Phong[3][3][10000], int a, int i, float u, float V[3])
{
    Phong[a][0][i] = (1 / u) * V[0];
    Phong[a][1][i] = (1 / u) * V[1];
    Phong[a][2][i] = (1 / u) * V[2];
}

void normalize4(float Phong[3][3][10000], int a, int i, float u)
{
    Phong[a][0][i] = (1 / u) * Phong[a][0][i];
    Phong[a][1][i] = (1 / u) * Phong[a][1][i];
    Phong[a][2][i] = (1 / u) * Phong[a][2][i];
}

float dot_product(float P1X, float P1Y, float P1Z, float P2X, float P2Y, float P2Z)
{
    float dot = (P1X * P2X) + (P1Y * P2Y) + (P1Z * P2Z);
    return dot;
}

void drawLinha(int V[3], float P_1[2], float P_2[2], float P_3[2], float P1x, float P1y, float P2x, float P2y)
{
    float auxX1 = P1x;
    float auxY1 = P1y;
    float auxX2 = P2x;
    float auxY2 = P2y;
    int dX = P2x - P1x;
    int dY = P2y - P1y;
    int eixoMaior = abs(dX) >= abs(dY) ? abs(dX) : abs(dY);
    int eixoMenor = abs(dX) < abs(dY) ? abs(dX) : abs(dY);
    int incrX0 = dX > 0 ? 1 : dX < 0 ? -1
                                     : 0;
    int incrY0 = dY > 0 ? 1 : dY < 0 ? -1
                                     : 0;
    int incrX1 = dX > 0 ? 1 : dX < 0 ? -1
                                     : 0;
    int incrY1 = dY > 0 ? 1 : dY < 0 ? -1
                                     : 0;

    if (abs(dX) < abs(dY))
    {
        eixoMaior = abs(dY);
        eixoMenor = abs(dX);
        incrX1 = 0;
        incrY1 = dY > 0 ? 1 : dY < 0 ? -1
                                     : 0;
    }
    int numerador = eixoMaior / 2;
    for (int i = 0; i <= eixoMaior; i++)
    {
        float ponto[2] = {auxX1, auxY1};
        float originalPoint[2][3];
        pOriginal(V, ponto, P_1, P_2, P_3, originalPoint);
        if (originalPoint[0][2] < zbuffer[(int)auxX1][(int)auxY1] && auxX1 >= 0 && auxX1 <= 800 && auxY1 >= 0 && auxY1 <= 450)
        {
            zbuffer[(int)auxX1][(int)auxY1] = (int)originalPoint[0][2];
            float vetorV[3] = {(originalPoint[0][0]) * (-1), (originalPoint[0][1]) * (-1), (originalPoint[0][2]) * (-1)};
            float normalV = sqrt(pow(vetorV[0], 2) + pow(vetorV[1], 2) + pow(vetorV[2], 2));
            normalize(vetorV, normalV);
            float vetorL[3] = {camera[0][13] - originalPoint[0][0], camera[1][13] - originalPoint[0][1], camera[2][13] - originalPoint[0][2]};
            float normalL = sqrt(pow(vetorL[0], 2) + pow(vetorL[1], 2) + pow(vetorL[2], 2));
            normalize(vetorL, normalL);
            float NL = dot_product(originalPoint[1][0], originalPoint[1][1], originalPoint[1][2], vetorL[0], vetorL[1], vetorL[2]);
            float NV = dot_product(originalPoint[1][0], originalPoint[1][1], originalPoint[1][2], vetorV[0], vetorV[1], vetorV[2]);
            float vetorR[3];
            for (int i = 0; i < 3; i++)
            {
                vetorR[i] = NL * 2 * originalPoint[1][i] - vetorL[i];
            }
            float RV = dot_product(vetorR[0], vetorR[1], vetorR[2], vetorV[0], vetorV[1], vetorV[2]);
            float Ia[3];
            float Id[3];
            float Is[3];
            for (int i = 0; i < 3; i++)
            {
                Ia[i] = camera[i][6] * camera[0][8];
                Id[i] = NL * camera[i][11] * camera[i][12] * camera[i][7];
                Is[i] = pow(RV, camera[0][10]) * camera[0][9] * camera[i][7];
            }

            if (NL < 0)
            {
                if (NV < 0)
                {
                    for (int i = 0; i < 3; i++)
                    {
                        originalPoint[1][i] *= -1;
                    }
                    NL = dot_product(originalPoint[1][0], originalPoint[1][1], originalPoint[1][2], vetorL[0], vetorL[1], vetorL[2]);
                    for (int i = 0; i < 3; i++)
                    {
                        vetorR[i] = (NL * 2 * originalPoint[1][i]) - vetorL[i];
                    }
                    RV = dot_product(vetorR[0], vetorR[1], vetorR[2], vetorV[0], vetorV[1], vetorV[2]);
                    for (int i = 0; i < 3; i++)
                    {
                        Id[i] = NL * camera[i][11] * camera[i][12] * camera[i][7];
                        Is[i] = pow(RV, camera[0][10]) * camera[0][9] * camera[i][7];
                    }
                }
                else
                {
                    Id[0] = 0;
                    Id[1] = 0;
                    Id[2] = 0;

                    Is[0] = 0;
                    Is[1] = 0;
                    Is[2] = 0;
                }
            }
            if (RV < 0)
            {
                Is[0] = 0;
                Is[1] = 0;
                Is[2] = 0;
            }
            float I[3] = {Ia[0] + Id[0] + Is[0], Ia[1] + Id[1] + Is[0], Ia[2] + Id[2] + Is[0]};
            if (I[0] > 255)
            {
                I[0] = 255;
            }
            if (I[1] > 255)
            {
                I[1] = 255;
            }
            if (I[2] > 255)
            {
                I[2] = 255;
            }
            DrawPixel(auxX1, auxY1, (Color){I[0], I[1], I[2], 255});
        }
        numerador += eixoMenor;
        if (numerador > eixoMaior)
        {
            numerador -= eixoMaior;
            auxX1 += incrX0;
            auxY1 += incrY0;
        }
        else
        {
            auxX1 += incrX1;
            auxY1 += incrY1;
        }
    }
}

void swap(float *a, float *b)
{
    float temp = *a;
    *a = *b;
    *b = temp;
}

void organizar(float pontoCima[2], float pontoMedio[2], float pontoBaixo[2])
{
    float ponto1X = pontoCima[0];
    float ponto1Y = pontoCima[1];
    float ponto2X = pontoMedio[0];
    float ponto2Y = pontoMedio[1];
    float ponto3X = pontoBaixo[0];
    float ponto3Y = pontoBaixo[1];

    if (ponto1Y > ponto3Y)
    {
        swap(&ponto3Y, &ponto1Y);
        swap(&ponto3X, &ponto1X);
    }
    if (ponto1Y > ponto2Y)
    {
        swap(&ponto2Y, &ponto1Y);
        swap(&ponto2X, &ponto1X);
    }
    if (ponto2Y > ponto3Y)
    {
        swap(&ponto3Y, &ponto2Y);
        swap(&ponto3X, &ponto2X);
    }

    if (ponto1Y == ponto3Y)
    {
        if (ponto1X > ponto3X)
        {
            swap(&ponto3Y, &ponto1Y);
            swap(&ponto3X, &ponto1X);
        }
    }
    if (ponto1Y == ponto2Y)
    {
        if (ponto1X > ponto2X)
        {
            swap(&ponto2Y, &ponto1Y);
            swap(&ponto2X, &ponto1X);
        }
    }
    if (ponto2Y == ponto3Y)
    {
        if (ponto2X > ponto3X)
        {
            swap(&ponto3Y, &ponto2Y);
            swap(&ponto3X, &ponto2X);
        }
    }

    pontoCima[0] = ponto1X;
    pontoCima[1] = ponto1Y;
    pontoMedio[0] = ponto2X;
    pontoMedio[1] = ponto2Y;
    pontoBaixo[0] = ponto3X;
    pontoBaixo[1] = ponto3Y;
}

void rasterizar(int V[3], float Pixel1x, float Pixel1y, float Pixel2x, float Pixel2y, float Pixel3x, float Pixel3y)
{
    float pontoAlto[2] = {Pixel1x, Pixel1y};
    float pontoIntermediario[2] = {Pixel2x, Pixel2y};
    float pontoMenor[2] = {Pixel3x, Pixel3y};
    organizar(pontoAlto, pontoIntermediario, pontoMenor);
    float Yanalisador = pontoAlto[1];
    float Xminimo = pontoAlto[0];
    float Xmaximo = pontoAlto[0];
    float Esquerda, Direita;
    float temp = (pontoIntermediario[0] - pontoAlto[0]);
    if (temp == 0)
    {
        Esquerda = 0;
    }
    else
    {
        Esquerda = (int)(pontoIntermediario[1] - pontoAlto[1]) / (pontoIntermediario[0] - pontoAlto[0]);
        if (Esquerda != 0)
        {
            Esquerda = (1 / Esquerda);
        }
    }
    temp = (pontoMenor[0] - pontoAlto[0]);
    if (temp == 0)
    {
        Direita = 0;
    }
    else
    {
        Direita = (int)(pontoMenor[1] - pontoAlto[1]) / (pontoMenor[0] - pontoAlto[0]);
        if (!Direita == 0)
        {
            Direita = (1 / Direita);
        }
    }
    if (pontoIntermediario[0] > pontoMenor[0])
    {
        float troca = Esquerda;
        Esquerda = Direita;
        Direita = troca;
    }
    if (Esquerda == 0 && Direita == 0)
    {
        Xmaximo = pontoMenor[0];
    }
    while (Yanalisador <= pontoIntermediario[1])
    {
        float Xminimo_ = (Xminimo + Esquerda);
        float Xmaximo_ = (Xmaximo + Direita);
        drawLinha(V, pontoAlto, pontoIntermediario, pontoMenor, round(Xminimo), Yanalisador, round(Xmaximo), Yanalisador);
        Yanalisador++;
        Xminimo = Xminimo_;
        Xmaximo = Xmaximo_;
    }
    float pontoD[2];
    int xDiff = pontoIntermediario[0] > pontoMenor[0] ? Xminimo - Esquerda : Xmaximo - Direita;
    pontoD[0] = xDiff;
    pontoD[1] = pontoIntermediario[1];
    Yanalisador = pontoMenor[1];
    Xminimo = Xmaximo = pontoMenor[0];
    temp = pontoIntermediario[0] - pontoMenor[0];
    if (temp == 0)
    {
        Esquerda = 0;
    }
    else
    {
        Esquerda = (pontoIntermediario[1] - pontoMenor[1]) / temp;
        if (Esquerda != 0)
        {
            Esquerda = 1 / Esquerda;
        }
    }
    temp = pontoD[0] - pontoMenor[0];
    if (temp == 0)
    {
        Direita = 0;
    }
    else
    {
        Direita = (pontoD[1] - pontoMenor[1]) / temp;
        if (Direita != 0)
        {
            Direita = 1 / Direita;
        }
    }
    if (pontoIntermediario[0] > pontoD[0])
    {
        float troca = Esquerda;
        Esquerda = Direita;
        Direita = troca;
    }
    if (Esquerda == 0 && Direita == 0)
    {
        Xmaximo = pontoIntermediario[0];
    }
    while (Yanalisador >= pontoIntermediario[1])
    {
        float Xminimo_ = Xminimo - Esquerda;
        float Xmaximo_ = Xmaximo - Direita;
        drawLinha(V, pontoAlto, pontoIntermediario, pontoMenor, round(Xminimo), Yanalisador, round(Xmaximo), Yanalisador);
        Yanalisador--;
        Xminimo = Xminimo_;
        Xmaximo = Xmaximo_;
    }
}

void desenhar(float teste[3][3][10000])
{
    for (int i = 0; i < teste[0][1][0]; i++)
    {
        int V[3] = {teste[2][0][i] - 1, teste[2][1][i] - 1, teste[2][2][i] - 1};
        rasterizar(V, teste[1][0][V[0]], teste[1][1][V[0]], teste[1][0][V[1]], teste[1][1][V[1]], teste[1][0][V[2]], teste[1][1][V[2]]);
    }
}

void coordenadaTela(float parametros[3][3][10000], float camera[3][14])
{
    int numVertices = parametros[0][0][0];
    float vertices[3][10000];
    for (int i = 0; i < numVertices; i++)
    {
        vertices[0][i] = parametros[1][0][i];
        vertices[1][i] = parametros[1][1][i];
        vertices[2][i] = parametros[1][2][i];
    }
    float C[3] = {camera[0][5], camera[1][5], camera[2][5]};
    float N[3] = {camera[0][0], camera[1][0], camera[2][0]};
    float V[3] = {camera[0][1], camera[1][1], camera[2][1]};
    float d = camera[0][2];
    float hx = camera[0][3];
    float hy = camera[0][4];
    for (int i = 0; i < numVertices; i++)
    {
        float P[3] = {vertices[0][i], vertices[1][i], vertices[2][i]};
        float calcV = (V[0] * N[0] + V[1] * N[1] + V[2] * N[2]) / (N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
        float Vortogonal[3] = {V[0] - (calcV * N[0]), V[1] - (calcV * N[1]), V[2] - (calcV * N[2])};
        float U[3] = {(N[1] * Vortogonal[2]) - (N[2] * Vortogonal[1]), (Vortogonal[0] * N[2]) - (Vortogonal[2] * N[0]), (N[0] * Vortogonal[1]) - (N[1] * Vortogonal[0])};
        float U_[3];
        float _U_ = sqrt(pow(U[0], 2) + pow(U[1], 2) + pow(U[2], 2));
        normalize2(U_, _U_, U);
        float V_[3];
        float _V_ = sqrt(pow(Vortogonal[0], 2) + pow(Vortogonal[1], 2) + pow(Vortogonal[2], 2));
        normalize2(V_, _V_, Vortogonal);
        float N_[3];
        float _N_ = sqrt(pow(N[0], 2) + pow(N[1], 2) + pow(N[2], 2));
        normalize2(N_, _N_, N);
        float P_[3] = {P[0] - C[0], P[1] - C[1], P[2] - C[2]};
        float Pvista[3];
        Pvista[0] = U_[0] * P_[0] + U_[1] * P_[1] + U_[2] * P_[2];
        Pvista[1] = V_[0] * P_[0] + V_[1] * P_[1] + V_[2] * P_[2];
        Pvista[2] = N_[0] * P_[0] + N_[1] * P_[1] + N_[2] * P_[2];
        for (int j = 0; j < 3; j++)
        {
            Phong[0][j][i] = Pvista[j];
        }
        float projPersp[2] = {d * (Pvista[0] / Pvista[2]), d * (Pvista[1] / Pvista[2])};
        float cordNormalizada[2] = {projPersp[0] / hx, projPersp[1] / hy};
        float pTela[2];
        pTela[0] = round(((cordNormalizada[0] + 1) / 2) * screenWidth + 0.5);
        pTela[1] = round(screenHeight - (((cordNormalizada[1] + 1) / 2) * screenHeight + 0.5));
        parametros[1][0][i] = pTela[0];
        parametros[1][1][i] = pTela[1];
    }
}

void valPhong(float teste[3][3][10000])
{
    int numVertices = teste[0][0][0];
    int numTriangulos = teste[0][1][0];
    float triangulos[3][10000];

    for (int i = 0; i < numTriangulos; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            triangulos[j][i] = teste[2][j][i];
        }

        int ver1 = triangulos[0][i] - 1;
        int ver2 = triangulos[1][i] - 1;
        int ver3 = triangulos[2][i] - 1;

        float vet1[3], vet2[3], normal[3];
        for (int j = 0; j < 3; j++)
        {
            vet1[j] = Phong[0][j][ver2] - Phong[0][j][ver1];
            vet2[j] = Phong[0][j][ver3] - Phong[0][j][ver1];
            normal[j] = (j == 0 ? vet1[1] * vet2[2] - vet1[2] * vet2[1]
                                : (j == 1 ? vet2[0] * vet1[2] - vet2[2] * vet1[0]
                                          : vet1[0] * vet2[1] - vet1[1] * vet2[0]));
        }

        float _normal_ = sqrt(pow(normal[0], 2) + pow(normal[1], 2) + pow(normal[2], 2));
        normalize3(Phong, 1, i, _normal_, normal);
    }

    for (int i = 0; i < numVertices; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Phong[2][j][i] = 0.0;
        }

        for (int j = 0; j < numTriangulos; j++)
        {
            if (triangulos[0][j] == i + 1 || triangulos[1][j] == i + 1 || triangulos[2][j] == i + 1)
            {
                for (int k = 0; k < 3; k++)
                {
                    Phong[2][k][i] += Phong[1][k][j];
                }
            }
        }

        float _normal_ = sqrt(pow(Phong[2][0][i], 2) + pow(Phong[2][1][i], 2) + pow(Phong[2][2][i], 2));
        normalize4(Phong, 2, i, _normal_);
    }
}

void reiniciarPhong()
{
    memset(Phong, 0, sizeof(Phong));
}

void reiniciarBuffer()
{
    memset(zbuffer, 999999, sizeof(zbuffer));
}

void dividirDados(char *nome, float dados[])
{
    int contador = 0;

    char *token = strtok(nome, " ");
    while (token != NULL && contador < 6)
    {
        dados[contador] = atof(token);
        contador++;
        token = strtok(NULL, " ");
    }
}

void Malha(float teste[3][3][10000])
{
    float dados[3];
    int contador = 0, contador1 = 0, contador2 = 0;

    FILE *fp = fopen(filename, "r");
    if (fp == NULL)
    {
        printf("Error: nao foi possivel ler o arquivo %s", filename);
        return;
    }

    char line[256];
    while (fgets(line, sizeof(line), fp))
    {
        sscanf(line, "%f %f %f", &dados[0], &dados[1], &dados[2]);

        if (contador == 0)
        {
            teste[0][0][0] = dados[0];
            teste[0][1][0] = dados[1];
        }
        else if (contador > 0 && contador <= teste[0][0][0])
        {
            teste[1][0][contador1] = dados[0];
            teste[1][1][contador1] = dados[1];
            teste[1][2][contador1] = dados[2];

            contador1++;
        }
        else if (contador > teste[0][0][0])
        {
            teste[2][0][contador2] = dados[0];
            teste[2][1][contador2] = dados[1];
            teste[2][2][contador2] = dados[2];

            contador2++;
        }

        contador++;
    }

    fclose(fp);
}

void visaoCamera(float teste[3][14])
{
    char filename[100] = "parametros/camera.txt";
    FILE *fp = fopen(filename, "r");
    if (fp == NULL)
    {
        printf("Error: cannot open file %s", filename);
        return;
    }
    int contador = 0;
    float *dados = malloc(10 * sizeof(float));

    const unsigned MAX_LENGTH = 256;
    char buffer[MAX_LENGTH];

    while (fgets(buffer, MAX_LENGTH, fp))
    {
        dividirDados(buffer, dados);

        teste[0][contador] = dados[2];
        teste[1][contador] = dados[3];
        teste[2][contador] = dados[4];

        contador++;
    }
    fclose(fp);
    free(dados);
}

int main(void)
{
    InitWindow(screenWidth, screenHeight, "Projeto Computação Gráfica");
    Malha(obj);
    visaoCamera(camera);
    reiniciarBuffer(zbuffer);
    coordenadaTela(obj, camera);
    valPhong(obj);

    while (!WindowShouldClose())
    {

        if (IsKeyPressed('R'))
        {
            Malha(obj);
            visaoCamera(camera);
            reiniciarBuffer();
            reiniciarPhong();
            coordenadaTela(obj, camera);
            valPhong(obj);
        }
        BeginDrawing();
        ClearBackground(BLACK);
        desenhar(obj);
        reiniciarBuffer();
        EndDrawing();
    }

    CloseWindow();
    return 0;
}