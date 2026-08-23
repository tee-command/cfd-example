#include <math.h>
#include <stdio.h>

// =============================================================================
// グローバル定数の宣言
// =============================================================================

// 計算領域の無次元長さ
const double LX = 2.0;
const double LY = 1.0;

// 格子点数
const int NX = 64;
const int NY = 32;

// 計算領域の開始座標と終点座標
const int X_MIN = 0;
const int Y_MIN = 0;
const int X_MAX = NX - 1;
const int Y_MAX = NY - 1;

// 格子間隔
const double DX = LX / NX;
const double DY = LY / NY;

// 物理量
const double RE = 1000;

// シミュレーション時間の設定
const double DT = 0.001;
const double SIMULATION_MAX_TIME = 10.0;

// SOR法の設定
const double SOR_OMEGA = 1.5; // 緩和係数（1〜2の間で調整、収束性を見て設定）
const int SOR_MAX_ITER = 1000;
const double SOR_EPS = 1e-8; // 収束判定閾値

// =============================================================================
// サブルーチン：出力
// =============================================================================

// 速度ノルムがしきい値以上なら画面に@や0で表示
void print_velocities(double u[NX][NY], double v[NX][NY]) {
  for (int y = Y_MAX; y >= Y_MIN; --y) {
    for (int x = X_MIN; x <= X_MAX; ++x) {
      const double velocity_norm = sqrt(u[x][y] * u[x][y] + v[x][y] * v[x][y]);
      // 速度ノルムが0.1以上なら@で表示
      if (velocity_norm >= 0.1) {
        printf("@");
        continue;
      }

      // 速度ノルムが0.05以上なら0で表示
      if (velocity_norm >= 0.05) {
        printf("0");
        continue;
      }

      // 速度ノルムが0.01以上ならoで表示
      if (velocity_norm >= 0.01) {
        printf("o");
        continue;
      }

      // 速度ノルムが0.01未満なら.で表示
      printf(".");
    }
    // 次の行へ
    printf("\n");
  }
  // 出力が終われば区切りを引く
  for (int i = X_MIN; i <= X_MAX; ++i) {
    printf("=");
  }
  printf("\n");
}

// =============================================================================
// サブルーチン：数学操作
// =============================================================================

// a (mod b)
int mod(int a, int b) { return (a + b) % b; }

// =============================================================================
// サブルーチン：初期化
// =============================================================================

// 配列u,v,pの要素を0埋めする
void zero_fill_velocities_and_pressure_fields(double u[NX][NY],
                                              double v[NX][NY],
                                              double p[NX][NY]) {
  for (int y = Y_MIN; y <= Y_MAX; ++y) {
    for (int x = X_MIN; x <= X_MAX; ++x) {
      u[x][y] = 0;
      v[x][y] = 0;
      p[x][y] = 0;
    }
  }
};

// 速度場を任意に設定
void init_velocities(double u[NX][NY], double v[NX][NY]) {
  u[NX / 2 - 4][NY / 2 - 2] = 1.5;
  u[NX / 2 - 4][NY / 2 - 1] = 1.5;
  u[NX / 2 - 4][NY / 2] = 1.5;
  u[NX / 2 - 4][NY / 2 + 1] = 1.5;
  u[NX / 2 - 4][NY / 2 + 2] = 1.5;

  u[NX / 2 - 3][NY / 2 - 2] = 1.5;
  u[NX / 2 - 3][NY / 2 - 1] = 1.5;
  u[NX / 2 - 3][NY / 2] = 1.5;
  u[NX / 2 - 3][NY / 2 + 1] = 1.5;
  u[NX / 2 - 3][NY / 2 + 2] = 1.5;

  u[NX / 2 - 2][NY / 2 - 2] = 1.5;
  u[NX / 2 - 2][NY / 2 - 1] = 1.5;
  u[NX / 2 - 2][NY / 2] = 1.5;
  u[NX / 2 - 2][NY / 2 + 1] = 1.5;
  u[NX / 2 - 2][NY / 2 + 2] = 1.5;

  u[NX / 2 - 1][NY / 2 - 2] = 1.5;
  u[NX / 2 - 1][NY / 2 - 1] = 1.5;
  u[NX / 2 - 1][NY / 2] = 1.5;
  u[NX / 2 - 1][NY / 2 + 1] = 1.5;
  u[NX / 2 - 1][NY / 2 + 2] = 1.5;

  u[NX / 2][NY / 2 - 2] = 1.5;
  u[NX / 2][NY / 2 - 1] = 1.5;
  u[NX / 2][NY / 2] = 1.5;
  u[NX / 2][NY / 2 + 1] = 1.5;
  u[NX / 2][NY / 2 + 2] = 1.5;

  u[NX / 2 + 1][NY / 2 - 2] = 1.5;
  u[NX / 2 + 1][NY / 2 - 1] = 1.5;
  u[NX / 2 + 1][NY / 2] = 1.5;
  u[NX / 2 + 1][NY / 2 + 1] = 1.5;
  u[NX / 2 + 1][NY / 2 + 2] = 1.5;

  u[NX / 2 + 2][NY / 2 - 2] = 1.5;
  u[NX / 2 + 2][NY / 2 - 1] = 1.5;
  u[NX / 2 + 2][NY / 2] = 1.5;
  u[NX / 2 + 2][NY / 2 + 1] = 1.5;
  u[NX / 2 + 2][NY / 2 + 2] = 1.5;

  u[NX / 2 + 3][NY / 2 - 2] = 1.5;
  u[NX / 2 + 3][NY / 2 - 1] = 1.5;
  u[NX / 2 + 3][NY / 2] = 1.5;
  u[NX / 2 + 3][NY / 2 + 1] = 1.5;
  u[NX / 2 + 3][NY / 2 + 2] = 1.5;

  u[NX / 2 + 3][NY / 2 - 2] = 1.5;
  u[NX / 2 + 3][NY / 2 - 1] = 1.5;
  u[NX / 2 + 3][NY / 2] = 1.5;
  u[NX / 2 + 3][NY / 2 + 1] = 1.5;
  u[NX / 2 + 3][NY / 2 + 2] = 1.5;
}

// =============================================================================
// サブルーチン：KKスキーム
// =============================================================================

// u.d(phi)/d(xi)の離散化
double kkscheme_x(const double u, const double phi_imm, const double phi_im,
                  const double phi_i, const double phi_ip,
                  const double phi_ipp) {
  return u * (-phi_ipp + 8 * phi_ip - 8 * phi_im + phi_imm) / (12 * DX) +
         fabs(u) * (phi_ipp - 4 * phi_ip + 6 * phi_i - 4 * phi_im + phi_imm) /
             (4 * DX);
}

double kkscheme_y(const double u, const double phi_imm, const double phi_im,
                  const double phi_i, const double phi_ip,
                  const double phi_ipp) {
  return u * (-phi_ipp + 8 * phi_ip - 8 * phi_im + phi_imm) / (12 * DY) +
         fabs(u) * (phi_ipp - 4 * phi_ip + 6 * phi_i - 4 * phi_im + phi_imm) /
             (4 * DY);
}

// =============================================================================
// サブルーチン：ラプラシアン
// =============================================================================

// 中心点の係数: -30/(12h^2) = -2.5/h^2
const double LAPLACIAN4_DIAG_X = -2.5 / (DX * DX);
const double LAPLACIAN4_DIAG_Y = -2.5 / (DY * DY);
const double LAPLACIAN4_DIAG = LAPLACIAN4_DIAG_X + LAPLACIAN4_DIAG_Y;

// 4次中心差分ラプラシアンの近傍項（中心点 phi[x][y] を除く）
double laplacian4_off_diag(const double phi[NX][NY], int x, int y) {
  const double dx2 = (-phi[mod(x - 2, NX)][y] + 16 * phi[mod(x - 1, NX)][y] +
                      16 * phi[mod(x + 1, NX)][y] - phi[mod(x + 2, NX)][y]) /
                     (12 * DX * DX);
  const double dy2 = (-phi[x][mod(y - 2, NY)] + 16 * phi[x][mod(y - 1, NY)] +
                      16 * phi[x][mod(y + 1, NY)] - phi[x][mod(y + 2, NY)]) /
                     (12 * DY * DY);
  return dx2 + dy2;
}

// フルのラプラシアン（RHSの陽的評価用）
double laplacian4(const double phi[NX][NY], int x, int y) {
  return laplacian4_off_diag(phi, x, y) + LAPLACIAN4_DIAG * phi[x][y];
}

// pの4次精度ラプラシアンの近傍項（中心点 p[x][y] を除く）
double laplacian4_off_diag_p(const double p[NX][NY], int x, int y) {
  const double dx2 = (-p[mod(x - 2, NX)][y] + 16 * p[mod(x - 1, NX)][y] +
                      16 * p[mod(x + 1, NX)][y] - p[mod(x + 2, NX)][y]) /
                     (12 * DX * DX);
  const double dy2 = (-p[x][mod(y - 2, NY)] + 16 * p[x][mod(y - 1, NY)] +
                      16 * p[x][mod(y + 1, NY)] - p[x][mod(y + 2, NY)]) /
                     (12 * DY * DY);
  return dx2 + dy2;
}

// =============================================================================
// サブルーチン：Navier-Stokes
// =============================================================================

void correct_velocities_by_pressure(double u[NX][NY], double v[NX][NY],
                                    const double u_star_x[NX][NY],
                                    const double u_star_y[NX][NY],
                                    const double p[NX][NY]) {
  for (int y = Y_MIN; y <= Y_MAX; ++y) {
    for (int x = X_MIN; x <= X_MAX; ++x) {
      // u[x][y] は p[x-1][y] と p[x][y] に挟まれた位置にある
      const double dpdx = (p[x][y] - p[mod(x - 1, NX)][y]) / DX;
      u[x][y] = u_star_x[x][y] - DT * dpdx;

      // v[x][y] は p[x][y-1] と p[x][y] に挟まれた位置にある
      const double dpdy = (p[x][y] - p[x][mod(y - 1, NY)]) / DY;
      v[x][y] = u_star_y[x][y] - DT * dpdy;
    }
  }
}

void solve_p(double p[NX][NY], const double u_star_x[NX][NY],
             const double u_star_y[NX][NY]) {
  // RHS = (RHO/DT) * div(u_star)   （スタガード格子の1階中心差分）
  double rhs[NX][NY];
  for (int y = Y_MIN; y <= Y_MAX; ++y) {
    for (int x = X_MIN; x <= X_MAX; ++x) {
      const double div = (u_star_x[mod(x + 1, NX)][y] - u_star_x[x][y]) / DX +
                         (u_star_y[x][mod(y + 1, NY)] - u_star_y[x][y]) / DY;
      rhs[x][y] = DT * div;
    }
  }

  // SOR法で Lap(p) = rhs を解く
  for (int iter = 0; iter < SOR_MAX_ITER; ++iter) {
    double max_residual = 0.0;

    for (int y = Y_MAX; y >= Y_MIN; --y) {
      for (int x = X_MAX; x >= X_MIN; --x) {
        const double p_gs =
            (rhs[x][y] - laplacian4_off_diag_p(p, x, y)) / LAPLACIAN4_DIAG;
        const double res = p_gs - p[x][y];
        p[x][y] += SOR_OMEGA * res;
        max_residual = fmax(max_residual, fabs(res));
      }
    }

    if (max_residual < SOR_EPS)
      break;
  }
}

void solve_u_and_v_star(double u_star_x[NX][NY], double u_star_y[NX][NY],
                        const double u[NX][NY], const double v[NX][NY],
                        const double advection_x[NX][NY],
                        const double advection_y[NX][NY],
                        const double advection_x_old[NX][NY],
                        const double advection_y_old[NX][NY]) {
  // RHS = u^n + dt*(1.5*advection^n - 0.5*advection^{n-1})
  //       + 0.5*dt*Lap(u^n)/RE　（すべて既知量）
  double rhs_x[NX][NY];
  double rhs_y[NX][NY];
  for (int y = Y_MIN; y <= Y_MAX; ++y) {
    for (int x = X_MIN; x <= X_MAX; ++x) {
      rhs_x[x][y] =
          u[x][y] +
          DT * (1.5 * advection_x[x][y] - 0.5 * advection_x_old[x][y]) +
          0.5 * DT * laplacian4(u, x, y) / RE;
      rhs_y[x][y] =
          v[x][y] +
          DT * (1.5 * advection_y[x][y] - 0.5 * advection_y_old[x][y]) +
          0.5 * DT * laplacian4(v, x, y) / RE;
    }
  }

  for (int y = Y_MIN; y <= Y_MAX; ++y) {
    for (int x = X_MIN; x <= X_MAX; ++x) {
      u_star_x[x][y] = u[x][y];
      u_star_y[x][y] = v[x][y];
    }
  }

  const double A = 0.5 * DT / RE;
  const double denom = 1.0 - A * LAPLACIAN4_DIAG;

  for (int iter = 0; iter < SOR_MAX_ITER; ++iter) {
    double max_residual = 0.0;
    for (int y = Y_MAX; y >= Y_MIN; --y) {
      for (int x = X_MAX; x >= X_MIN; --x) {
        const double u_gs =
            (rhs_x[x][y] + A * laplacian4_off_diag(u_star_x, x, y)) / denom;
        const double res_u = u_gs - u_star_x[x][y];
        u_star_x[x][y] += SOR_OMEGA * res_u;
        max_residual = fmax(max_residual, fabs(res_u));

        const double v_gs =
            (rhs_y[x][y] + A * laplacian4_off_diag(u_star_y, x, y)) / denom;
        const double res_v = v_gs - u_star_y[x][y];
        u_star_y[x][y] += SOR_OMEGA * res_v;
        max_residual = fmax(max_residual, fabs(res_v));
      }
    }
    if (max_residual < SOR_EPS)
      break;
  }
}

void solve_advection(double advection_x[NX][NY], double advection_y[NX][NY],
                     const double u[NX][NY], const double v[NX][NY]) {
  for (int y = Y_MAX; y >= Y_MIN; --y) {
    for (int x = X_MAX; x >= X_MIN; --x) {
      // uの定義点でのududx, vdudy
      const double ududx_on_u =
          kkscheme_x(u[x][y], u[mod(x - 2, NX)][y], u[mod(x - 1, NX)][y],
                     u[x][y], u[mod(x + 1, NX)][y], u[mod(x + 2, NX)][y]);
      const double vdudy_on_u =
          kkscheme_y((v[x][y] + v[mod(x - 1, NX)][y] + v[x][mod(y + 1, NY)] +
                      v[mod(x - 1, NX)][mod(y + 1, NY)]) /
                         4,
                     u[x][mod(y - 2, NY)], u[x][mod(y - 1, NY)], u[x][y],
                     u[x][mod(y + 1, NY)], u[x][mod(y + 2, NY)]);

      // vの定義点でのudvdx, vdvdy
      const double udvdx_on_v =
          kkscheme_x((u[x][y] + u[mod(x + 1, NX)][y] + u[x][mod(y - 1, NY)] +
                      u[mod(x + 1, NX)][mod(y - 1, NY)]) /
                         4,
                     v[mod(x - 2, NX)][y], v[mod(x - 1, NX)][y], v[x][y],
                     v[mod(x + 1, NX)][y], v[mod(x + 2, NX)][y]);
      const double vdvdy_on_v =
          kkscheme_y(v[x][y], v[x][mod(y - 2, NY)], v[x][mod(y - 1, NY)],
                     v[x][y], v[x][mod(y + 1, NY)], v[x][mod(y + 2, NY)]);

      advection_x[x][y] = -ududx_on_u - vdudy_on_u;
      advection_y[x][y] = -udvdx_on_v - vdvdy_on_v;
    }
  }
}

void update_fields(double u[NX][NY], double v[NX][NY], double p[NX][NY],
                   double advection_x_old[NX][NY],
                   double advection_y_old[NX][NY], const int step) {
  // 時刻nの移流項
  double advection_x[NX][NY];
  double advection_y[NX][NY];
  solve_advection(advection_x, advection_y, u, v);

  // 初回は１次精度Euler陽解法として移流項を求める
  if (step == 0) {
    for (int y = Y_MIN; y <= Y_MAX; ++y) {
      for (int x = X_MIN; x <= X_MAX; ++x) {
        advection_x_old[x][y] = advection_x[x][y];
        advection_y_old[x][y] = advection_y[x][y];
      }
    }
  }

  double u_star_x[NX][NY];
  double u_star_y[NX][NY];
  solve_u_and_v_star(u_star_x, u_star_y, u, v, advection_x, advection_y,
                     advection_x_old, advection_y_old);

  solve_p(p, u_star_x, u_star_y);
  correct_velocities_by_pressure(u, v, u_star_x, u_star_y, p);

  // 次ステップのAB2のために今回の移流項を保存
  for (int y = Y_MIN; y <= Y_MAX; ++y) {
    for (int x = X_MIN; x <= X_MAX; ++x) {
      advection_x_old[x][y] = advection_x[x][y];
      advection_y_old[x][y] = advection_y[x][y];
    }
  }
}

// =============================================================================
// メイン関数
//
// 計算領域　：x-y周期境界条件
// 速度場　　：スタッガード
// 圧力場　　：コロケート
// 時間進行　：フラクショナルステップ
// 離散化　　：４次精度中心差分
// 陰解法　　：SOR法
// =============================================================================
int main() {
  // 速度場・圧力場を管理する配列の宣言
  double u[NX][NY];
  double v[NX][NY];
  double p[NX][NY];

  // AB2用：1ステップ前の移流項を保持する配列
  double advection_x_old[NX][NY];
  double advection_y_old[NX][NY];
  int step = 0;

  // 速度圧力の初期化
  zero_fill_velocities_and_pressure_fields(u, v, p);
  init_velocities(u, v);

  while (1) {
    // 目標時間に達したらループ終了
    const double elapsed_time = step * DT;
    if (elapsed_time > SIMULATION_MAX_TIME) {
      break;
    }

    // 1000 step ごとに画面出力
    if (step % 1000 == 0) {
      print_velocities(u, v);
    }

    // 時間進行
    update_fields(u, v, p, advection_x_old, advection_y_old, step);
    step++;
  }

  // 画面出力
  print_velocities(u, v);

  return 0;
}
