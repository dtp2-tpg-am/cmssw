#ifndef ANALTYPDEFS_H
#define ANALTYPDEFS_H
#include "constants.h"
#include <stdint.h>

/* Posibles calidades de la trayectoria:
   NOPATH => No es una trayectoria válida
   LOWQGHOST => Igual que LOWQ con múltiples casos simultáneos
   LOWQ   => Es una potencial trayectoria pero sólo formada por 3 puntos
   HIGHQGHOST => Igual que HIGHQ con múltiples casos simultáneos
   HIGHQ  => Es una trayectoria válida con 4 puntos alineados (4 celdas)
OB*/
typedef enum {NOPATH = 0, LOWQGHOST, LOWQ, HIGHQGHOST, HIGHQ} MP_QUALITY;
// Tipos de lateralidad de traza de partícula al pasar por una celda
typedef enum {LEFT=0, RIGHT, NONE} LATERAL_CASES;
typedef struct {
  uint32_t rawId;
  float t0;
  float x;
  float tanPhi;
  float phi;
  float phiB;
  float chi2;
  int quality;
  int wi1;
  int tdc1;
  int wi2;
  int tdc2;
  int wi3;
  int tdc3;
  int wi4;
  int tdc4;
  int wi5;
  int tdc5;
  int wi6;
  int tdc6;
  int wi7;
  int tdc7;
  int wi8;
  int tdc8;
} metaPrimitive;
#endif
