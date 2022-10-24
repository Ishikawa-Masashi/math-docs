---
title: 回転行列、クォータニオン(四元数)、オイラー角の相互変換
tags: 3D 数学
author: aa_debdeb
slide: false
---
3D CGでは回転を表現するのに回転行列、クォータニオン、オイラー角の3つがよく使われます。この記事では、その3つの表現を相互に変換する方法について解説します。

# 準備

回転行列 $a \ne 0$ $ \boldsymbol{R}$の各要素を次のように表します。

$$ \displaystyle
\boldsymbol{R} =
\left(
  \begin{array}{ccc}
    m_{00} & m_{01} & m_{02} \\
    m_{10} & m_{11} & m_{12} \\
    m_{20} & m_{21} & m_{22}
  \end{array}
\right)\\
$$

クォータニオン$\boldsymbol{q}$を次のように表します。

$$
\boldsymbol{q} =
q_w + q_xi + q_yj + q_zk
=
\left[
  \begin{array}{c}
    q_x \\
    q_y \\
    q_z \\
    q_w
  \end{array}
\right] \\
i^2 + j^2 + k ^2 = -1\\
ij = k\\
jk = i\\
ki = j
$$

オイラー角のXYZ各軸に対する回転角度をそれぞれ$\theta_x$、$\theta_y$、$\theta_z$と表現します。オイラー角は各軸に対する回転を適用する順番により最終的な結果が変わります。この記事では回転順を後から適用されるものから順に軸を並べたXYZというような形式で表現します。例えばXYZの場合、Z軸、Y軸、X軸という順に回転を適用していきます。


# オイラー角から回転行列

オイラー角から回転行列への変換です。

X軸に対する回転行列$\boldsymbol{R_x}$、Y軸に対する回転行列$\boldsymbol{R_y}$、Z軸に対する回転行列$\boldsymbol{R_z}$はそれぞれ次のようになります。

```math
\begin{eqnarray}
\boldsymbol{R_x} &=&
\left(
  \begin{array}{ccc}
    1 & 0 & 0 \\
    0 & \cos\theta_{x} & -\sin\theta_{x} \\
    0 & \sin\theta_{x} & \cos\theta_{x}
  \end{array}
\right) \\
\boldsymbol{R_y} &=&
\left(
  \begin{array}{ccc}
    \cos\theta_{y} & 0 & \sin\theta_{y} \\
    0 & 1 & 0 \\
    -\sin\theta_{y}  & 0 & \cos\theta_{y}
  \end{array}
\right) \\
\boldsymbol{R_z} &=&
\left(
  \begin{array}{ccc}
    \cos\theta_{z} & -\sin\theta_{z} & 0 \\
    \sin\theta_{z} & \cos\theta_{z} & 0 \\
    0 & 0 & 1
  \end{array}
\right)
\end{eqnarray}
```

この3つの行列$\boldsymbol{R_x}$、$\boldsymbol{R_y}$、$\boldsymbol{R_z}$の積を求めることで、オイラー角を回転行列に変換することができます。先述したようにオイラー角では回転順が最終結果に影響するので、各回転順での変換を求めていきます。

### 回転順XYZ

```math
\boldsymbol{R_{xyz}} = \boldsymbol{R_x}\boldsymbol{R_y}\boldsymbol{R_z} = 
\left(
  \begin{array}{ccc}
    \cos\theta_{y}\cos\theta_{z} & -\cos\theta_{y}\sin\theta_{z} & \sin\theta_{y} \\
    \sin\theta_{x}\sin\theta_{y}\cos\theta_{z} + \cos\theta_{x}\sin\theta_{z} & -\sin\theta_{x}\sin\theta_{y}\sin\theta_{z} + \cos\theta_{x}\cos\theta_{z} & -\sin\theta_{x}\cos\theta_{y} \\
    -\cos\theta_{x}\sin\theta_{y}\cos\theta_{z} + \sin\theta_{x}\sin\theta_{z} & \cos\theta_{x}\sin\theta_{y}\sin\theta_{z} + \sin\theta_{x}\cos\theta_{z} & \cos\theta_{x}\cos\theta_{y}
  \end{array}
\right)
```

### 回転順XZY

```math
\boldsymbol{R_{xzy}} = \boldsymbol{R_x}\boldsymbol{R_z}\boldsymbol{R_y} = 
\left(
  \begin{array}{ccc}
    \cos\theta_{y}\cos\theta_{z} & -\sin\theta_{z} & \sin\theta_{y}\cos\theta_{z} \\
    \cos\theta_{x}\cos\theta_{y}\sin\theta_{z} + \sin\theta_{x}\sin\theta_{y} & \cos\theta_{x}\cos\theta_{z} & \cos\theta_{x}\sin\theta_{y}\sin\theta_{z} - \sin\theta_{x}\cos\theta_{y} \\
    \sin\theta_{x}\cos\theta_{y}\sin\theta_{z} - \cos\theta_{x}\sin\theta_{y} & \sin\theta_{x}\cos\theta_{z} & \sin\theta_{x}\sin\theta_{y}\sin\theta_{z} + \cos\theta_{x}\cos\theta_{y}
  \end{array}
\right)\\
```

### 回転順YXZ

```math
\boldsymbol{R_{yxz}} = \boldsymbol{R_y}\boldsymbol{R_x}\boldsymbol{R_z} = 
\left(
  \begin{array}{ccc}
    \sin\theta_{x}\sin\theta_{y}\sin\theta_{z} + \cos\theta_{y}\cos\theta_{z} & \sin\theta_{x}\sin\theta_{y}\cos\theta_{z} -\cos\theta_{y}\sin\theta_{z} & \cos\theta_{x}\sin\theta_{y} \\
    \cos\theta_{x}\sin\theta_{z} & \cos\theta_{x}\cos\theta_{z} & -\sin\theta_{x} \\
    \sin\theta_{x}\cos\theta_{y}\sin\theta_{z} -\sin\theta_{y}\cos\theta_{z} & \sin\theta_{x}\cos\theta _{y}\cos\theta_{z} + \sin\theta_{y}\sin\theta_{z} & \cos\theta_{x}\cos\theta_{y}
  \end{array}
\right)\\
```

### 回転順YZX

```math
\boldsymbol{R_{yzx}} = \boldsymbol{R_y}\boldsymbol{R_z}\boldsymbol{R_x} = 
\left(
  \begin{array}{ccc}
    \cos\theta_{y}\cos\theta_{z} & -\cos\theta_{x}\cos\theta_{y}\sin\theta_{z} + \sin\theta_{x}\sin\theta_{y} & \sin\theta_{x}\cos\theta_{y}\sin\theta_{z} + \cos\theta_{x}\sin\theta_{y} \\
   \sin\theta_{z} & \cos\theta_{x}\cos\theta_{z} & -\sin\theta_{x}\cos\theta_{z} \\
    -\sin\theta_{y}\cos\theta_{z} & \cos\theta_{x}\sin\theta_{y}\sin\theta_{z} + \sin\theta_{x}\cos\theta_{y} & -\sin\theta_{x}\sin\theta_{y}\sin\theta_{z} + \cos\theta_{x}\cos\theta_{y}
  \end{array}
\right)
```

### 回転順ZXY

```math
\boldsymbol{R_{zxy}} = \boldsymbol{R_z}\boldsymbol{R_x}\boldsymbol{R_y} = 
\left(
  \begin{array}{ccc}
    -\sin\theta_{x}\sin\theta_{y}\sin\theta_{z} + \cos\theta_{y}\cos\theta_{z} & -\cos\theta_{x}\sin\theta_{z} & \sin\theta_{x}\cos\theta_{y}\sin\theta_{z} + \sin\theta_{y}\cos\theta_{z} \\
    \sin\theta_{x}\sin\theta_{y}\cos\theta_{z} + \cos\theta_{y}\sin\theta_{z} & \cos\theta_{x}\cos\theta_{z} & -\sin\theta_{x}\cos\theta_{y}\cos\theta_{z} + \sin\theta_{y}\sin\theta_{z} \\
    -\cos\theta_{x}\sin\theta_{y} & \sin\theta_{x} & \cos\theta_{x}\cos\theta_{y}
  \end{array}
\right)
```

### 回転順ZYX

```math
\boldsymbol{R_{zyx}} = \boldsymbol{R_z}\boldsymbol{R_y}\boldsymbol{R_x} = 
\left(
  \begin{array}{ccc}
    \cos\theta_{y}\cos\theta_{z} & \sin\theta_{x}\sin\theta_{y}\cos\theta_{z} - \cos\theta_{x}\sin\theta_{z} & \cos\theta_{x}\sin\theta_{y}\cos\theta_{z} + \sin\theta_{x}\sin\theta_{z} \\
    \cos\theta_{y}\sin\theta_{z} & \sin\theta_{x}\sin\theta_{y}\sin\theta_{z} + \cos\theta_{x}\cos\theta_{z} & \cos\theta_{x}\sin\theta_{y}\sin\theta_{z} -\sin\theta_{x}\cos\theta_{z} \\
    -\sin\theta_{y} & \sin\theta_{x}\cos\theta _{y} & \cos\theta_{x}\cos\theta_{y}
  \end{array}
\right)\\
```

# 回転行列からオイラー角

回転行列からオイラー角への変換です。

先ほど求めたオイラー角から回転行列への変換をもとにして、回転行列からオイラー角への変換を求めます。

### 回転順XYZ

回転順XYZのオイラー角から回転行列への変換は次のようになります。

```math
\boldsymbol{R_{xyz}} = 
\left(
  \begin{array}{ccc}
    \cos\theta_{y}\cos\theta_{z} & -\cos\theta_{y}\sin\theta_{z} & \sin\theta_{y} \\
    \sin\theta_{x}\sin\theta_{y}\cos\theta_{z} + \cos\theta_{x}\sin\theta_{z} & -\sin\theta_{x}\sin\theta_{y}\sin\theta_{z} + \cos\theta_{x}\cos\theta_{z} & -\sin\theta_{x}\cos\theta_{y} \\
    -\cos\theta_{x}\sin\theta_{y}\cos\theta_{z} + \sin\theta_{x}\sin\theta_{z} & \cos\theta_{x}\sin\theta_{y}\sin\theta_{z} + \sin\theta_{x}\cos\theta_{z} & \cos\theta_{x}\cos\theta_{y}
  \end{array}
\right)
```

ここから、$\theta_y$が次のように求まります。

```math
m_{02} =  \sin\theta_{y} \\
\theta_{y} = \arcsin(m_{02})
```

また、$\cos\theta_y \neq 0$のとき、$\theta_x$、$\theta_z$が次のように求まります。

```math
\begin{eqnarray}
m_{12} &=& -\sin\theta_x\cos\theta_y \\
\sin\theta_x &=& -\frac{m_{12}}{\cos\theta_y} \\
m_{22} &=& \cos\theta_x\cos\theta_y \\
\cos\theta_x &=& \frac{m_{22}}{\cos\theta_y} \\
\tan\theta_x &=& \frac{\sin\theta_x}{\cos\theta_x} = -\frac{m_{12}}{m_{22}} \\
\theta_x &=& \arctan(-\frac{m_{12}}{m_{22}})
\end{eqnarray}
```

```math
\begin{eqnarray}
m_{00} &=& \cos\theta_y\cos\theta_z \\
\cos\theta_z &=& \frac{m_{00}}{\cos\theta_y}  \\
m_{01} &=& -\cos\theta_y\sin\theta_z \\
\sin\theta_z &=& -\frac{m_{01}}{\cos\theta_y}  \\
\tan\theta_z &=& \frac{\sin\theta_z}{\cos\theta_z} = -\frac{m_{01}}{m_{00}} \\
\theta_z &=& \arctan(-\frac{m_{01}}{m_{00}})
\end{eqnarray}
```

$\cos\theta_y = 0$のとき、ジンバルロックが発生してX軸とZ軸の回転が同軸での回転となります。$\theta_z=0$と仮定すると、回転行列$\boldsymbol{R_{xyz}}$は次のようになります。

```math
\boldsymbol{R_{xyz}} = 
\left(
  \begin{array}{ccc}
    0 & 0 & \sin\theta_{y} \\
    \sin\theta_{x}\sin\theta_{y} & \cos\theta_{x} & 0 \\
    -\cos\theta_{x}\sin\theta_{y} & \sin\theta_{x} & 0
  \end{array}
\right)
```

ここから$\theta_x$が次のように求まります。

```math
\tan\theta_x = \frac{\sin\theta_x}{\cos\theta_x} = \frac{m_{21}}{m_{11}} \\
\theta_x = \arctan(\frac{m_{21}}{m_{11}})
```

結果として、回転順XYZのときの回転行列からオイラー角への変換は次のようになります。

```math
\begin{eqnarray}
\theta_x &=& \left\{ \begin{array}{ll}
  \arctan(-\frac{m_{12}}{m_{22}}) & (cos\theta_y \neq 0) \\
  \arctan(\frac{m_{21}}{m_{11}}) & (otherwise)
\end{array} \right. \\
\theta_y &=& \arcsin(m_{02}) \\
\theta_z &=& \left\{ \begin{array}{ll}
  \arctan(-\frac{m_{01}}{m_{00}}) & (cos\theta_y \neq 0) \\
  0 & (otherwise)
\end{array} \right. \\
\end{eqnarray}
```

以下、各回転順での回転行列からオイラー角への変換を同様の方法で求めていきます。

### 回転順XZY

```math
\boldsymbol{R_{xzy}} = 
\left(
  \begin{array}{ccc}
    \cos\theta_{y}\cos\theta_{z} & -\sin\theta_{z} & \sin\theta_{y}\cos\theta_{z} \\
    \cos\theta_{x}\cos\theta_{y}\sin\theta_{z} + \sin\theta_{x}\sin\theta_{y} & \cos\theta_{x}\cos\theta_{z} & \cos\theta_{x}\sin\theta_{y}\sin\theta_{z} - \sin\theta_{x}\cos\theta_{y} \\
    \sin\theta_{x}\cos\theta_{y}\sin\theta_{z} - \cos\theta_{x}\sin\theta_{y} & \sin\theta_{x}\cos\theta_{z} & \sin\theta_{x}\sin\theta_{y}\sin\theta_{z} + \cos\theta_{x}\cos\theta_{y}
  \end{array}
\right)\\
```

より、

```math
\begin{eqnarray}
m_{01} &=& -\sin\theta_z \\
\theta_z &=& \arcsin(-m_{01})
\end{eqnarray}
```

$\cos\theta_z \neq 0$のとき、

```math
\begin{eqnarray}
m_{21} &=& \sin\theta_x\cos\theta_z \\
\sin\theta_x &=& \frac{m_{21}}{\cos\theta_z} \\
m_{11} &=& \cos\theta_x\cos\theta_z \\
\cos\theta_x &=& \frac{m_{11}}{\cos\theta_z} \\
\tan\theta_x &=& \frac{\sin\theta_x}{\cos\theta_x} = \frac{m_{21}}{m_{11}} \\
\theta_x &=& \arctan(\frac{m_{21}}{m_{11}}) \\
m_{02} &=& \sin\theta_y\cos\theta_z \\
\sin\theta_y &=& \frac{m_{02}}{\cos\theta_z}  \\
m_{00} &=& \cos\theta_y\cos\theta_z \\
\cos\theta_y &=& \frac{m_{00}}{\cos\theta_z}  \\
\tan\theta_y &=& \frac{\sin\theta_y}{\cos\theta_y} = \frac{m_{02}}{m_{00}} \\
\theta_y &=& \arctan(\frac{m_{02}}{m_{00}})
\end{eqnarray}
```

$\cos\theta_z = 0$のとき、ジンバルロックにより$\theta_y = 0$と仮定すると、

```math
\boldsymbol{R_{xzy}} = 
\left(
  \begin{array}{ccc}
    0 & -\sin\theta_z & 0 \\
   \cos\theta_x\sin\theta_z & 0 & -\sin\theta_x \\
   \sin\theta_x\sin\theta_z & 0 & \cos\theta_x
  \end{array}
\right)
```

よって、

```math
\begin{eqnarray}
\tan\theta_x &=& \frac{\sin\theta_x}{\cos\theta_x} = -\frac{m_{12}}{m_{22}} \\
\theta_x &=& \arctan(-\frac{m_{12}}{m_{22}})
\end{eqnarray}
```

結果、

```math
\begin{eqnarray}
\theta_x &=& \left\{ \begin{array}{ll}
  \arctan(\frac{m_{21}}{m_{11}}) & (cos\theta_z \neq 0) \\
  \arctan(-\frac{m_{12}}{m_{22}}) & (otherwise)
\end{array} \right. \\
\theta_y &=& \left\{ \begin{array}{ll}
  \arctan(\frac{m_{02}}{m_{00}}) & (cos\theta_z \neq 0) \\
  0 & (otherwise)
\end{array} \right. \\
\theta_z &=& \arcsin(-m_{01}) \\
\end{eqnarray}
```

### 回転順YXZ

```math
\boldsymbol{R_{yxz}} =
\left(
  \begin{array}{ccc}
    \sin\theta_{x}\sin\theta_{y}\sin\theta_{z} + \cos\theta_{y}\cos\theta_{z} & \sin\theta_{x}\sin\theta_{y}\cos\theta_{z} -\cos\theta_{y}\sin\theta_{z} & \cos\theta_{x}\sin\theta_{y} \\
    \cos\theta_{x}\sin\theta_{z} & \cos\theta_{x}\cos\theta_{z} & -\sin\theta_{x} \\
    \sin\theta_{x}\cos\theta_{y}\sin\theta_{z} -\sin\theta_{y}\cos\theta_{z} & \sin\theta_{x}\cos\theta _{y}\cos\theta_{z} + \sin\theta_{y}\sin\theta_{z} & \cos\theta_{x}\cos\theta_{y}
  \end{array}
\right)\\
```

より、

```math
\begin{eqnarray}
m_{12} &=& -\sin\theta_x \\
\theta_x &=& \arcsin(-m_{12})
\end{eqnarray}
```

$\cos\theta_x \neq 0$のとき、

```math
\begin{eqnarray}
m_{02} &=& \cos\theta_x\sin\theta_y \\
\sin\theta_y &=& \frac{m_{02}}{\cos\theta_x} \\
m_{22} &=& \cos\theta_x\cos\theta_y \\
\cos\theta_y &=& \frac{m_{22}}{\cos\theta_x} \\
\tan\theta_y &=& \frac{\sin\theta_y}{\cos\theta_y} = \frac{m_{02}}{m_{22}} \\
\theta_y &=& \arctan(\frac{m_{02}}{m_{22}}) \\
m_{10} &=& \cos\theta_x\sin\theta_z \\
\sin\theta_z &=& \frac{m_{10}}{\cos\theta_x}  \\
m_{11} &=& \cos\theta_x\cos\theta_z \\
\cos\theta_z &=& \frac{m_{11}}{\cos\theta_x}  \\
\tan\theta_z &=& \frac{\sin\theta_z}{\cos\theta_z} = \frac{m_{10}}{m_{11}} \\
\theta_z &=& \arctan(\frac{m_{10}}{m_{11}})
\end{eqnarray}
```

$\cos\theta_x = 0$のとき、ジンバルロックにより$\theta_z = 0$と仮定すると、

```math
\boldsymbol{R_{yxz}} = 
\left(
  \begin{array}{ccc}
    \cos\theta_y & \sin\theta_x\sin\theta_y & 0 \\
    0 & 0 & -\sin\theta_x \\
    -\sin\theta_y & \sin\theta_x\sin\theta_y & 0
  \end{array}
\right)
```

よって、

```math
\begin{eqnarray}
\tan\theta_y &=& \frac{\sin\theta_y}{\cos\theta_y} = -\frac{m_{20}}{m_{00}} \\
\theta_y &=& \arctan(-\frac{m_{20}}{m_{00}})
\end{eqnarray}
```

結果、

```math
\begin{eqnarray}
\theta_x &=& \arcsin(-m_{12}) \\
\theta_y &=& \left\{ \begin{array}{ll}
  \arctan(\frac{m_{02}}{m_{22}}) & (cos\theta_x \neq 0) \\
  \arctan(-\frac{m_{20}}{m_{00}}) & (otherwise)
\end{array} \right. \\
\theta_z &=& \left\{ \begin{array}{ll}
  \arctan(\frac{m_{10}}{m_{11}}) & (cos\theta_x \neq 0) \\
  0 & (otherwise)
\end{array} \right. \\
\end{eqnarray}
```

### 回転順YZX

```math
\boldsymbol{R_{yzx}} = 
\left(
  \begin{array}{ccc}
    \cos\theta_{y}\cos\theta_{z} & -\cos\theta_{x}\cos\theta_{y}\sin\theta_{z} + \sin\theta_{x}\sin\theta_{y} & \sin\theta_{x}\cos\theta_{y}\sin\theta_{z} + \cos\theta_{x}\sin\theta_{y} \\
   \sin\theta_{z} & \cos\theta_{x}\cos\theta_{z} & -\sin\theta_{x}\cos\theta_{z} \\
    -\sin\theta_{y}\cos\theta_{z} & \cos\theta_{x}\sin\theta_{y}\sin\theta_{z} + \sin\theta_{x}\cos\theta_{y} & -\sin\theta_{x}\sin\theta_{y}\sin\theta_{z} + \cos\theta_{x}\cos\theta_{y}
  \end{array}
\right)
```

より、

```math
\begin{eqnarray}
m_{10} &=& \sin\theta_z \\
\theta_z &=& \arcsin(m_{10})
\end{eqnarray}
```

$\cos\theta_z \neq 0$のとき、

```math
\begin{eqnarray}
m_{12} &=& -\sin\theta_x\cos\theta_z \\
\sin\theta_x &=& -\frac{m_{12}}{\cos\theta_z} \\
m_{11} &=& \cos\theta_x\cos\theta_z \\
\cos\theta_x &=& \frac{m_{11}}{\cos\theta_z} \\
\tan\theta_x &=& \frac{\sin\theta_x}{\cos\theta_x} = -\frac{m_{12}}{m_{11}} \\
\theta_x &=& \arctan(-\frac{m_{12}}{m_{11}}) \\
m_{20} &=& -\sin\theta_y\cos\theta_z \\
\sin\theta_y &=& -\frac{m_{20}}{\cos\theta_z}  \\
m_{00} &=& \cos\theta_y\cos\theta_z \\
\cos\theta_y &=& \frac{m_{00}}{\cos\theta_z}  \\
\tan\theta_y &=& \frac{\sin\theta_y}{\cos\theta_y} = -\frac{m_{20}}{m_{00}} \\
\theta_y &=& \arctan(-\frac{m_{20}}{m_{00}})
\end{eqnarray}
```

$\cos\theta_z = 0$のとき、ジンバルロックにより$\theta_x = 0$と仮定すると、

```math
\boldsymbol{R_{yzx}} = 
\left(
  \begin{array}{ccc}
    0 & -\cos\theta_y\sin\theta_z & \sin\theta_y \\
    \sin\theta_z & 0 & 0 \\
    0 & \sin\theta_y\sin\theta_z & \cos\theta_y
  \end{array}
\right)
```

よって、

```math
\begin{eqnarray}
\tan\theta_y &=& \frac{\sin\theta_y}{\cos\theta_y} = \frac{m_{02}}{m_{22}} \\
\theta_y &=& \arctan(\frac{m_{02}}{m_{22}})
\end{eqnarray}
```

結果、

```math
\begin{eqnarray}
\theta_x &=& \left\{ \begin{array}{ll}
  \arctan(-\frac{m_{12}}{m_{11}}) & (cos\theta_z \neq 0) \\
  0 & (otherwise)
\end{array} \right. \\
\theta_y &=& \left\{ \begin{array}{ll}
  \arctan(-\frac{m_{20}}{m_{00}}) & (cos\theta_z \neq 0) \\
  \arctan(\frac{m_{02}}{m_{22}}) & (otherwise)
\end{array} \right. \\
\theta_z &=& \arcsin(m_{10}) \\
\end{eqnarray}
```

### 回転順ZXY

```math
\boldsymbol{R_{zxy}} = 
\left(
  \begin{array}{ccc}
    -\sin\theta_{x}\sin\theta_{y}\sin\theta_{z} + \cos\theta_{y}\cos\theta_{z} & -\cos\theta_{x}\sin\theta_{z} & \sin\theta_{x}\cos\theta_{y}\sin\theta_{z} + \sin\theta_{y}\cos\theta_{z} \\
    \sin\theta_{x}\sin\theta_{y}\cos\theta_{z} + \cos\theta_{y}\sin\theta_{z} & \cos\theta_{x}\cos\theta_{z} & -\sin\theta_{x}\cos\theta_{y}\cos\theta_{z} + \sin\theta_{y}\sin\theta_{z} \\
    -\cos\theta_{x}\sin\theta_{y} & \sin\theta_{x} & \cos\theta_{x}\cos\theta_{y}
  \end{array}
\right)
```

より、

```math
\begin{eqnarray}
m_{21} &=& \sin\theta_x \\
\theta_x &=& \arcsin(m_{21})
\end{eqnarray}
```

$\cos\theta_x \neq 0$のとき、

```math
\begin{eqnarray}
m_{20} &=& -\cos\theta_x\sin\theta_y \\
\sin\theta_y &=& -\frac{m_{20}}{\cos\theta_x} \\
m_{22} &=& \sin\theta_x\cos\theta_y \\
\cos\theta_y &=& \frac{m_{22}}{\sin\theta_x} \\
\tan\theta_y &=& \frac{\sin\theta_y}{\cos\theta_y} = -\frac{m_{20}}{m_{22}} \\
\theta_y &=& \arctan(-\frac{m_{20}}{m_{21}}) \\
m_{01} &=& -\cos\theta_x\sin\theta_z \\
\sin\theta_z &=& -\frac{m_{01}}{\cos\theta_x}  \\
m_{11} &=& \cos\theta_x\cos\theta_z \\
\cos\theta_z &=& \frac{m_{11}}{\cos\theta_x}  \\
\tan\theta_z &=& \frac{\sin\theta_z}{\cos\theta_z} = -\frac{m_{01}}{m_{11}} \\
\theta_z &=& \arctan(-\frac{m_{01}}{m_{11}})
\end{eqnarray}
```

$\cos\theta_x = 0$のとき、ジンバルロックにより$\theta_y = 0$と仮定すると、

```math
\boldsymbol{R_{zxy}} = 
\left(
  \begin{array}{ccc}
    \cos\theta_z & 0 & \sin\theta_x\sin\theta_z \\
    \sin\theta_z & 0 & -\sin\theta_x\cos\theta_z \\
    0 & \sin\theta_x & 0
  \end{array}
\right)
```

よって、

```math
\begin{eqnarray}
\tan\theta_z &=& \frac{\sin\theta_z}{\cos\theta_z} = \frac{m_{10}}{m_{00}} \\
\theta_z &=& \arctan(\frac{m_{10}}{m_{00}})
\end{eqnarray}
```

結果、

```math
\begin{eqnarray}
\theta_x &=& \arcsin(m_{21}) \\
\theta_y &=& \left\{ \begin{array}{ll}
  \arctan(-\frac{m_{20}}{m_{22}}) & (cos\theta_x \neq 0) \\
  0 & (otherwise)
\end{array} \right. \\
\theta_z &=& \left\{ \begin{array}{ll}
  \arctan(-\frac{m_{01}}{m_{11}}) & (cos\theta_x \neq 0) \\
  \arctan(\frac{m_{10}}{m_{00}}) & (otherwise)
\end{array} \right. \\
\end{eqnarray}
```

### 回転順ZYX

```math
\boldsymbol{R_{zyx}} =
\left(
  \begin{array}{ccc}
    \cos\theta_{y}\cos\theta_{z} & \sin\theta_{x}\sin\theta_{y}\cos\theta_{z} - \cos\theta_{x}\sin\theta_{z} & \cos\theta_{x}\sin\theta_{y}\cos\theta_{z} + \sin\theta_{x}\sin\theta_{z} \\
    \cos\theta_{y}\sin\theta_{z} & \sin\theta_{x}\sin\theta_{y}\sin\theta_{z} + \cos\theta_{x}\cos\theta_{z} & \cos\theta_{x}\sin\theta_{y}\sin\theta_{z} -\sin\theta_{x}\cos\theta_{z} \\
    -\sin\theta_{y} & \sin\theta_{x}\cos\theta _{y} & \cos\theta_{x}\cos\theta_{y}
  \end{array}
\right)\\
```

より、

```math
\begin{eqnarray}
m_{20} &=& -\sin\theta_y \\
\theta_y &=& \arcsin(-m_{20})
\end{eqnarray}
```

$\cos\theta_y \neq 0$のとき、

```math
\begin{eqnarray}
m_{21} &=& \sin\theta_x\cos\theta_y \\
\sin\theta_x &=& \frac{m_{21}}{\cos\theta_y} \\
m_{22} &=& \cos\theta_x\cos\theta_y \\
\cos\theta_x &=& \frac{m_{22}}{\sin\theta_y} \\
\tan\theta_x &=& \frac{\sin\theta_x}{\cos\theta_x} = \frac{m_{21}}{m_{22}} \\
\theta_x &=& \arctan(\frac{m_{21}}{m_{22}}) \\
m_{10} &=& \cos\theta_y\sin\theta_z \\
\sin\theta_z &=& \frac{m_{10}}{\cos\theta_y}  \\
m_{00} &=& \cos\theta_y\cos\theta_z \\
\cos\theta_z &=& \frac{m_{00}}{\cos\theta_y}  \\
\tan\theta_z &=& \frac{\sin\theta_z}{\cos\theta_z} = \frac{m_{10}}{m_{00}} \\
\theta_z &=& \arctan(\frac{m_{10}}{m_{00}})
\end{eqnarray}
```

$\cos\theta_y = 0$のとき、ジンバルロックにより$\theta_x = 0$と仮定すると、

```math
\boldsymbol{R_{zyx}} = 
\left(
  \begin{array}{ccc}
    0 & -\sin\theta_z & \sin\theta_y\cos\theta_z \\
    0 & \cos\theta_z & \sin\theta_y\sin\theta_z \\
    -\sin\theta_y & 0 & 0
  \end{array}
\right)
```

よって、

```math
\begin{eqnarray}
\tan\theta_z &=& \frac{\sin\theta_z}{\cos\theta_z} = -\frac{m_{01}}{m_{11}} \\
\theta_z &=& \arctan(-\frac{m_{01}}{m_{11}})
\end{eqnarray}
```

結果、

```math
\begin{eqnarray}
\theta_x &=& \left\{ \begin{array}{ll}
  \arctan(\frac{m_{21}}{m_{22}}) & (cos\theta_y \neq 0) \\
  0 & (otherwise)
\end{array} \right. \\
\theta_y &=& \arcsin(-m_{20}) \\
\theta_z &=& \left\{ \begin{array}{ll}
  \arctan(\frac{m_{10}}{m_{00}}) & (cos\theta_y \neq 0) \\
  \arctan(-\frac{m_{01}}{m_{11}}) & (otherwise)
\end{array} \right. \\
\end{eqnarray}
```


# オイラー角からクォータニオン

オイラー角からクォータニオンへの変換です。

単位ベクトル$\vec{n}$を軸として$\theta$だけ回転させるクォータニオンは次のようになります。

```math
\vec{n} = (n_x, n_y, v_z), \|\vec{n}\| = 1\\
\boldsymbol{q} =
\left[
  \begin{array}{c}
    n_x\sin\frac{\theta}{2} \\
    n_y\sin\frac{\theta}{2} \\
    n_z\sin\frac{\theta}{2} \\
    \cos\frac{\theta}{2}
  \end{array}
\right]
```

そのため、X軸、Y軸、Z軸周りの回転を表すクォータニオン$\boldsymbol{q_x}$、$\boldsymbol{q_y}$、$\boldsymbol{q_z}$はそれぞれ次のようになります。

```math
\boldsymbol{q_x} =
\left[
  \begin{array}{c}
    \sin\frac{\theta_x}{2} \\
    0 \\
    0 \\
    \cos\frac{\theta_x}{2}
  \end{array}
\right],
\boldsymbol{q_y} =
\left[
  \begin{array}{c}
    0 \\
    \sin\frac{\theta_y}{2} \\
    0 \\
    \cos\frac{\theta_y}{2}
  \end{array}
\right],
\boldsymbol{q_z} =
\left[
  \begin{array}{c}
    0 \\
    0 \\
    \sin\frac{\theta_z}{2} \\
    \cos\frac{\theta_z}{2}
  \end{array}
\right]
```

この3つのクォータニオン$\boldsymbol{q_x}$、$\boldsymbol{q_y}$、$\boldsymbol{q_z}$の積を求めることで、オイラー角からクォータニオンへの変換を求めることができます。

### 回転順XYZ

```math
\boldsymbol{q_x}\boldsymbol{q_y}\boldsymbol{q_z} = 
\left[
  \begin{array}{c}
    \cos\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\sin\frac{\theta_z}{2} + \sin\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\cos\frac{\theta_z}{2} \\
    -\sin\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\sin\frac{\theta_z}{2} + \cos\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\cos\frac{\theta_z}{2} \\
    \cos\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\sin\frac{\theta_z}{2} + \sin\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\cos\frac{\theta_z}{2} \\
    -\sin\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\sin\frac{\theta_z}{2} + \cos\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\cos\frac{\theta_z}{2}
  \end{array}
\right]
```

### 回転順XZY

```math
\boldsymbol{q_x}\boldsymbol{q_z}\boldsymbol{q_y} = 
\left[
  \begin{array}{c}
    -\cos\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\sin\frac{\theta_z}{2} + \sin\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\cos\frac{\theta_z}{2} \\
    \cos\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\cos\frac{\theta_z}{2} -\sin\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\sin\frac{\theta_z}{2} \\
    \sin\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\cos\frac{\theta_z}{2} + \cos\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\sin\frac{\theta_z}{2} \\
    \sin\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\sin\frac{\theta_z}{2} + \cos\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\cos\frac{\theta_z}{2}
  \end{array}
\right]
```

### 回転順YXZ

```math
\boldsymbol{q_y}\boldsymbol{q_x}\boldsymbol{q_z} = 
\left[
  \begin{array}{c}
    \cos\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\sin\frac{\theta_z}{2} + \sin\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\cos\frac{\theta_z}{2} \\
    -\sin\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\sin\frac{\theta_z}{2} + \cos\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\cos\frac{\theta_z}{2} \\
    \cos\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\sin\frac{\theta_z}{2} - \sin\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\cos\frac{\theta_z}{2} \\
    \sin\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\sin\frac{\theta_z}{2} + \cos\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\cos\frac{\theta_z}{2}
  \end{array}
\right]
```

### 回転順YZX

```math
\boldsymbol{q_y}\boldsymbol{q_z}\boldsymbol{q_x} = 
\left[
  \begin{array}{c}
    \sin\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\cos\frac{\theta_z}{2} + \cos\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\sin\frac{\theta_z}{2} \\
    \sin\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\sin\frac{\theta_z}{2} + \cos\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\cos\frac{\theta_z}{2} \\
    -\sin\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\cos\frac{\theta_z}{2} + \cos\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\sin\frac{\theta_z}{2} \\
    -\sin\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\sin\frac{\theta_z}{2} + \cos\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\cos\frac{\theta_z}{2}
  \end{array}
\right]
```

### 回転順ZXY

```math
\boldsymbol{q_z}\boldsymbol{q_x}\boldsymbol{q_y} = 
\left[
  \begin{array}{c}
    -\cos\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\sin\frac{\theta_z}{2} + \sin\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\cos\frac{\theta_z}{2} \\
    \cos\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\cos\frac{\theta_z}{2} + \sin\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\sin\frac{\theta_z}{2} \\
    \sin\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\cos\frac{\theta_z}{2} + \cos\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\sin\frac{\theta_z}{2} \\
    -\sin\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\sin\frac{\theta_z}{2} + \cos\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\cos\frac{\theta_z}{2}
  \end{array}
\right]
```

### 回転順ZYX

```math
\boldsymbol{q_z}\boldsymbol{q_y}\boldsymbol{q_x} = 
\left[
  \begin{array}{c}
    \sin\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\cos\frac{\theta_z}{2} - \cos\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\sin\frac{\theta_z}{2} \\
    \sin\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\sin\frac{\theta_z}{2} + \cos\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\cos\frac{\theta_z}{2} \\
    -\sin\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\cos\frac{\theta_z}{2} + \cos\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\sin\frac{\theta_z}{2} \\
    \sin\frac{\theta_x}{2}\sin\frac{\theta_y}{2}\sin\frac{\theta_z}{2} + \cos\frac{\theta_x}{2}\cos\frac{\theta_y}{2}\cos\frac{\theta_z}{2}
  \end{array}
\right]
```

# クォータニオンから回転行列

クォータニオンから回転行列への変換です。

単位ベクトル$\vec{n}$を軸として$\theta$だけ回転させる回転行列は次のようになります。

```math
\vec{n} = (n_x, n_y, v_z), \|\vec{n}\| = 1\\
\boldsymbol{R} =
\left(
  \begin{array}{ccc}
    \cos\theta + n_x^2(1 - \cos\theta) & n_xn_y(1 - \cos\theta) - n_z\sin\theta & n_xn_z(1 - \cos\theta) + n_y\sin\theta \\
    n_xn_y(1 - \cos\theta) + n_z\sin\theta & \cos\theta + n_y^2(1 - \cos\theta) & n_yn_z(1 - \cos\theta) - n_x\sin\theta \\
    n_xn_z(1 - \cos\theta) - n_y\sin\theta & n_yn_z(1 - cos\theta) + n_x\sin\theta & \cos\theta + n_z^2(1 - \cos\theta) 
  \end{array}
\right)
```

また、オイラー角からクォータニオンへの変換で紹介したように、任意軸に対する回転を表すクォータニオンは次のようになります。

```math
\vec{n} = (n_x, n_y, v_z), \|\vec{n}\| = 1\\
\boldsymbol{q} =
\left[
  \begin{array}{c}
    q_x \\
    q_y \\
    q_z \\
    q_w
  \end{array}
\right]
=
\left[
  \begin{array}{c}
    n_x\sin\frac{\theta}{2} \\
    n_y\sin\frac{\theta}{2} \\
    n_z\sin\frac{\theta}{2} \\
    \cos\frac{\theta}{2}
  \end{array}
\right]
```


下に示した三角関数の2倍角の公式を用いることで、行列の各要素をクォータニオンの各要素を用いた値に置き換えることができます。

```math
\begin{eqnarray}
\sin\theta &=& 2\sin\frac{\theta}{2}\cos\frac{\theta}{2}\\
\cos\theta &=& 2\cos^2\frac{\theta}{2} - 1 = 1 - 2\sin^2\frac{\theta}{2}
\end{eqnarray}
```

行列の各要素の置き換えは次のようになります。

```math
\begin{eqnarray}
m_{00} &=& \cos\theta + n_x^2(1 - \cos\theta) \\
&=& (2\cos^2\frac{\theta}{2} - 1) + n_x^2\{1 - (1 - 2\sin^2\frac{\theta}{2})\} \\
&=& 2\cos^2\frac{\theta}{2} + 2n_x^2\sin^2\frac{\theta}{2} - 1 \\
&=& 2q_w^2 + 2q_x^2 - 1 \\
m_{01} &=& n_xn_y(1 - \cos\theta) - n_z\sin\theta \\
&=& n_xn_y\{1 - (1 - 2\sin^2\frac{\theta}{2})\} - n_z(2\sin\frac{\theta}{2}\cos\frac{\theta}{2}) \\
&=& 2n_xn_y\sin^2\frac{\theta}{2} - 2n_z\sin\frac{\theta}{2}\cos\frac{\theta}{2} \\
&=& 2q_xq_y - 2q_zq_w \\
m_{02} &=& n_xn_z(1 - \cos\theta) + n_y\sin\theta \\
&=& n_xn_z\{1 - (1 - 2\sin^2\frac{\theta}{2})\} + n_y(2\sin\frac{\theta}{2}\cos\frac{\theta}{2}) \\
&=& 2n_xn_z\sin^2\frac{\theta}{2} + 2n_y\sin\frac{\theta}{2}\cos\frac{\theta}{2} \\
&=& 2q_xq_z + 2q_yq_w \\
m_{10} &=& n_xn_y(1 - \cos\theta) + n_z\sin\theta \\
&=& n_xn_y\{1 - (1 - 2\sin^2\frac{\theta}{2})\} + n_z(2\sin\frac{\theta}{2}\cos\frac{\theta}{2}) \\
&=& 2n_xn_y\sin^2\frac{\theta}{2} + 2n_z\sin\frac{\theta}{2}\cos\frac{\theta}{2} \\
&=& 2q_xq_y + 2q_zq_w \\
m_{11} &=& \cos\theta + n_y^2(1 - \cos\theta) \\
&=& (2\cos^2\frac{\theta}{2} - 1) + n_y^2\{1 - (1 - 2\sin^2\frac{\theta}{2})\} \\
&=& 2\cos^2\frac{\theta}{2} + 2n_y^2\sin^2\frac{\theta}{2} - 1\\
&=& 2q_w^2 + 2q_y^2 - 1\\
m_{12} &=& n_yn_z(1 - \cos\theta) - n_x\sin\theta \\
&=& n_yn_z\{1 - (1 - 2\sin^2\frac{\theta}{2})\} - n_x(2\sin\frac{\theta}{2}\cos\frac{\theta}{2}) \\
&=& 2n_yn_z\sin^2\frac{\theta}{2} - 2n_x\sin\frac{\theta}{2}\cos\frac{\theta}{2} \\
&=& 2q_yq_z - 2q_xq_w \\
m_{20} &=& n_xn_z(1 - \cos\theta) - n_y\sin\theta \\
&=& n_xn_z\{1 - (1 - 2\sin^2\frac{\theta}{2})\} - n_y(2\sin\frac{\theta}{2}\cos\frac{\theta}{2}) \\
&=& 2n_xn_z\sin^2\frac{\theta}{2} - 2n_y\sin\frac{\theta}{2}\cos\frac{\theta}{2}\\
&=& 2q_xq_z - 2q_yq_w \\
m_{21} &=& n_yn_z(1 - cos\theta) + n_x\sin\theta \\
&=& n_yn_z\{1 - (1 - 2\sin^2\frac{\theta}{2})\} + n_x(2\sin\frac{\theta}{2}\cos\frac{\theta}{2}) \\
&=& 2n_yn_z\sin^2\frac{\theta}{2} + 2n_x\sin\frac{\theta}{2}\cos\frac{\theta}{2} \\
&=& 2q_yq_z + 2q_xq_w \\
m_{22} &=& \cos\theta + n_z^2(1 - \cos\theta) \\
&=& (2\cos^2\frac{\theta}{2} - 1) + n_z^2\{1 - (1 - 2\sin^2\frac{\theta}{2})\} \\
&=& 2\cos^2\frac{\theta}{2} + 2n_z^2\sin^2\frac{\theta}{2} - 1\\
&=& 2q_w^2 + 2q_z^2 - 1 
\end{eqnarray}
```

結果、クォータニオンから回転行列への変換は次のようになります。

```math
\boldsymbol{R} =
\left(
  \begin{array}{ccc}
    2q_w^2 + 2q_x^2 - 1 & 2q_xq_y - 2q_zq_w & 2q_xq_z + 2q_yq_w \\
    2q_xq_y + 2q_zq_w & 2q_w^2 + 2q_y^2 - 1 & 2q_yq_z - 2q_xq_w \\
    2q_xq_z - 2q_yq_w & 2q_yq_z + 2q_xq_w & 2q_w^2 + 2q_z^2 - 1
  \end{array}
\right)
```

# 回転行列からクォータニオン

回転行列からクォータニオンへの変換です。

さきほど、クォータニオンから回転行列への変換が次のようになることを確認しました。

```math
\boldsymbol{R} =
\left(
  \begin{array}{ccc}
    2q_w^2 + 2q_x^2 - 1 & 2q_xq_y - 2q_zq_w & 2q_xq_z + 2q_yq_w \\
    2q_xq_y + 2q_zq_w & 2q_w^2 + 2q_y^2 - 1 & 2q_yq_z - 2q_xq_w \\
    2q_xq_z - 2q_yq_w & 2q_yq_z + 2q_xq_w & 2q_w^2 + 2q_z^2 - 1
  \end{array}
\right)
```

各要素の和・差をとります。

```math
\begin{eqnarray}
m_{10} + m_{01} &=& (2q_xq_y + 2q_zq_w) + (2q_xq_y - 2q_zq_w) = 4q_xq_y \\
m_{10} - m_{01} &=& (2q_xq_y + 2q_zq_w) - (2q_xq_y - 2q_zq_w) = 4q_zq_w \\
m_{02} + m_{20} &=& (2q_xq_z + 2q_yq_w) + (2q_xq_z - 2q_yq_w) = 4q_xq_z \\
m_{02} - m_{20} &=& (2q_xq_z + 2q_yq_w) - (2q_xq_z - 2q_yq_w) = 4q_yq_w \\
m_{21} + m_{12} &=& (2q_yq_z + 2q_xq_w) + (2q_yq_z - 2q_xq_w) = 4q_yq_z \\
m_{21} - m_{12} &=& (2q_yq_z + 2q_xq_w) - (2q_yq_z - 2q_xq_w) = 4q_xq_w \\
\end{eqnarray}
```

これを整理すると次のようになり、$q_x$、$q_y$,、$q_z$、$q_w$のうち1つをもとに他の値も求まることがわかります。

```math
q_x = \frac{m_{10} + m_{01}}{4q_y} = \frac{m_{02} + m_{20}}{4q_z} = \frac{m_{21} - m_{12}}{4q_w} \\
q_y = \frac{m_{10} + m_{01}}{4q_x} = \frac{m_{21} + m_{12}}{4q_z} = \frac{m_{02} - m_{20}}{4q_w} \\
q_z = \frac{m_{02} + m_{20}}{4q_x} = \frac{m_{21} + m_{12}}{4q_y} = \frac{m_{10} - m_{01}}{4q_w} \\
q_w = \frac{m_{21} - m_{12}}{4q_x} = \frac{m_{02} - m_{20}}{4q_y} = \frac{m_{10} - m_{01}}{4q_z}
```

$q_x$、$q_y$,、$q_z$、$q_w$の値は次のように求めることができます。

```math
\begin{eqnarray}
m_{00} - m_{11} - m_{22} &=& (2q_w^2 + 2q_x^2 - 1) - (2q_w^2 + 2q_y^2 - 1) - (2q_w^2 + 2q_z^2 - 1) \\
q_x^2 &=&  \frac{m_{00} - m_{11} - m_{22} + 1}{4} \\
q_x  &=&  \pm\frac{\sqrt{m_{00} - m_{11} - m_{22} + 1}}{2}
\end{eqnarray}
```

```math
\begin{eqnarray}
-m_{00} + m_{11} - m_{22} &=& -(2q_w^2 + 2q_x^2 - 1) + (2q_w^2 + 2q_y^2 - 1) - (2q_w^2 + 2q_z^2 - 1) \\
q_y^2 &=&  \frac{-m_{00} + m_{11} - m_{22} + 1}{4} \\
q_y &=&  \pm\frac{\sqrt{-m_{00} + m_{11} - m_{22} + 1}}{2}
\end{eqnarray}
```

```math
\begin{eqnarray}
-m_{00} - m_{11} + m_{22} &=& -(2q_w^2 + 2q_x^2 - 1) - (2q_w^2 + 2q_y^2 - 1) + (2q_w^2 + 2q_z^2 - 1) \\
q_z^2 &=& \frac{-m_{00} - m_{11} + m_{22} + 1}{4} \\
q_z &=& \pm\frac{\sqrt{-m_{00} - m_{11} + m_{22} + 1}}{2}
\end{eqnarray}
```

```math
\begin{eqnarray}
m_{00} + m_{11} + m_{22} &=& (2q_w^2 + 2q_x^2 - 1) + (2q_w^2 + 2q_y^2 - 1) + (2q_w^2 + 2q_z^2 - 1) \\
q_w^2 &=& \frac{m_{00} + m_{11} + m_{22} + 1}{2} \\
q_w &=& \pm\frac{\sqrt{m_{00} + m_{11} + m_{22} + 1}}{2}
\end{eqnarray}
```

$q_x$、$q_y$,、$q_z$、$q_w$の各解をもとに、回転行列からクォータニオンへの変換を整理すると次のようになります。解が負の場合は省略していますが、負の場合は各解を符号反転したものとなり、同じ回転になります。

```math
\boldsymbol{q} =
\left[
  \begin{array}{c}
    \frac{\sqrt{m_{00} - m_{11} - m_{22} + 1}}{2} \\
    \frac{m_{10} + m_{01}}{4q_x} \\
    \frac{m_{02} + m_{20}}{4q_x} \\
    \frac{m_{21} - m_{12}}{4q_x}
  \end{array}
\right]
=
\left[
  \begin{array}{c}
    \frac{\sqrt{m_{00} - m_{11} - m_{22} + 1}}{2} \\
    \frac{m_{10} + m_{01}}{2\sqrt{m_{00} - m_{11} - m_{22} + 1}} \\
    \frac{m_{02} + m_{20}}{2\sqrt{m_{00} - m_{11} - m_{22} + 1}} \\
    \frac{m_{21} - m_{12}}{2\sqrt{m_{00} - m_{11} - m_{22} + 1}}
  \end{array}
\right]
```

```math
\boldsymbol{q} =
\left[
  \begin{array}{c}
    \frac{m_{10} + m_{01}}{4q_y} \\
    \frac{\sqrt{-m_{00} + m_{11} - m_{22} + 1}}{2} \\
    \frac{m_{21} + m_{12}}{4q_y} \\
    \frac{m_{02} - m_{20}}{4q_y}
  \end{array}
\right]
=
\left[
  \begin{array}{c}
     \frac{m_{10} + m_{01}}{2\sqrt{-m_{00} + m_{11} - m_{22} + 1}} \\
    \frac{\sqrt{-m_{00} + m_{11} - m_{22} + 1}}{2} \\
    \frac{m_{21} + m_{12}}{2\sqrt{-m_{00} + m_{11} - m_{22} + 1}} \\
    \frac{m_{02} - m_{20}}{2\sqrt{-m_{00} + m_{11} - m_{22} + 1}}
  \end{array}
\right]
```

```math
\boldsymbol{q} =
\left[
  \begin{array}{c}
    \frac{m_{02} + m_{20}}{4q_z} \\
    \frac{m_{21} + m_{12}}{4q_z} \\
    \frac{\sqrt{-m_{00} - m_{11} + m_{22} + 1}}{2} \\
    \frac{m_{10} - m_{01}}{4q_z}
  \end{array}
\right]
=
\left[
  \begin{array}{c}
    \frac{m_{02} + m_{20}}{2\sqrt{-m_{00} - m_{11} + m_{22} + 1}} \\
    \frac{m_{21} + m_{12}}{2\sqrt{-m_{00} - m_{11} + m_{22} + 1}} \\
    \frac{\sqrt{-m_{00} - m_{11} + m_{22} + 1}}{2} \\
    \frac{m_{10} - m_{01}}{\sqrt{-m_{00} - m_{11} + m_{22} + 1}}
  \end{array}
\right]
```

```math
\boldsymbol{q} =
\left[
  \begin{array}{c}
    \frac{m_{21} - m_{12}}{4q_w} \\
    \frac{m_{02} - m_{20}}{4q_w} \\
    \frac{m_{10} - m_{01}}{4q_w} \\
    \frac{\sqrt{m_{00} + m_{11} + m_{22} + 1}}{2}
  \end{array}
\right]
=
\left[
  \begin{array}{c}
    \frac{m_{21} - m_{12}}{2\sqrt{m_{00} + m_{11} + m_{22} + 1}} \\
    \frac{m_{02} - m_{20}}{2\sqrt{m_{00} + m_{11} + m_{22} + 1}} \\
    \frac{m_{10} - m_{01}}{2\sqrt{m_{00} + m_{11} + m_{22} + 1}} \\
    \frac{\sqrt{m_{00} + m_{11} + m_{22} + 1}}{2}
  \end{array}
\right]
```

このように回転行列からクォータニオンへの変換では複数の変換候補がありますが、実務的にはゼロ除算が発生しないなど安定して解が求まる候補を選択します。

# クォータニオンからオイラー角

クォータニオンからオイラー角への変換です。

ここでは、先述したクォータニオンから回転行列への変換と回転行列からオイラー角への変換を利用します。クォータニオンから回転行列への変換は次のようになります。

```math
\boldsymbol{R} =
\left(
  \begin{array}{ccc}
    2q_w^2 + 2q_x^2 - 1 & 2q_xq_y - 2q_zq_w & 2q_xq_z + 2q_yq_w \\
    2q_xq_y + 2q_zq_w & 2q_w^2 + 2q_y^2 - 1 & 2q_yq_z - 2q_xq_w \\
    2q_xq_z - 2q_yq_w & 2q_yq_z + 2q_xq_w & 2q_w^2 + 2q_z^2 - 1
  \end{array}
\right)
```

この各要素を回転行列からオイラー角への変換に適用することで、クォータニオンからオイラー角への変換を求めます。

### 回転順XYZ

```math
\begin{eqnarray}
\theta_x &=& \left\{ \begin{array}{ll}
  \arctan(-\frac{m_{12}}{m_{22}}) & (cos\theta_y \neq 0) \\
  \arctan(\frac{m_{21}}{m_{11}}) & (otherwise)
\end{array} \right. \\
\theta_y &=& \arcsin(m_{02}) \\
\theta_z &=& \left\{ \begin{array}{ll}
  \arctan(-\frac{m_{01}}{m_{00}}) & (cos\theta_y \neq 0) \\
  0 & (otherwise)
\end{array} \right. \\
\end{eqnarray}
```

より、

```math
\begin{eqnarray}
\theta_x &=& \left\{ \begin{array}{ll}
  \arctan(-\frac{2q_yq_z - 2q_xq_w}{2q_w^2 + 2q_z^2 - 1}) & (cos\theta_y \neq 0) \\
  \arctan(\frac{2q_yq_z + 2q_xq_w}{2q_w^2 + 2q_y^2 - 1}) & (otherwise)
\end{array} \right. \\
\theta_y &=& \arcsin(2q_xq_z + 2q_yq_w) \\
\theta_z &=& \left\{ \begin{array}{ll}
  \arctan(-\frac{2q_xq_y - 2q_zq_w}{2q_w^2 + 2q_x^2 - 1}) & (cos\theta_y \neq 0) \\
  0 & (otherwise)
\end{array} \right. \\
\end{eqnarray}
```

### 回転順XZY

```math
\begin{eqnarray}
\theta_x &=& \left\{ \begin{array}{ll}
  \arctan(\frac{m_{21}}{m_{11}}) & (cos\theta_z \neq 0) \\
  \arctan(-\frac{m_{12}}{m_{22}}) & (otherwise)
\end{array} \right. \\
\theta_y &=& \left\{ \begin{array}{ll}
  \arctan(\frac{m_{02}}{m_{00}}) & (cos\theta_z \neq 0) \\
  0 & (otherwise)
\end{array} \right. \\
\theta_z &=& \arcsin(-m_{01}) \\
\end{eqnarray}
```

より、

```math
\begin{eqnarray}
\theta_x &=& \left\{ \begin{array}{ll}
  \arctan(\frac{2q_yq_z + 2q_xq_w}{2q_w^2 + 2q_y^2 - 1}) & (cos\theta_z \neq 0) \\
  \arctan(-\frac{2q_yq_z - 2q_xq_w}{2q_w^2 + 2q_z^2 - 1}) & (otherwise)
\end{array} \right. \\
\theta_y &=& \left\{ \begin{array}{ll}
  \arctan(\frac{2q_xq_z + 2q_yq_w}{2q_w^2 + 2q_x^2 - 1}) & (cos\theta_z \neq 0) \\
  0 & (otherwise)
\end{array} \right. \\
\theta_z &=& \arcsin(-(2q_xq_y - 2q_zq_w)) \\
\end{eqnarray}
```

### 回転順YXZ

```math
\begin{eqnarray}
\theta_x &=& \arcsin(-m_{12}) \\
\theta_y &=& \left\{ \begin{array}{ll}
  \arctan(\frac{m_{02}}{m_{22}}) & (cos\theta_x \neq 0) \\
  \arctan(-\frac{m_{20}}{m_{00}}) & (otherwise)
\end{array} \right. \\
\theta_z &=& \left\{ \begin{array}{ll}
  \arctan(\frac{m_{10}}{m_{11}}) & (cos\theta_x \neq 0) \\
  0 & (otherwise)
\end{array} \right. \\
\end{eqnarray}
```

より、

```math
\begin{eqnarray}
\theta_x &=& \arcsin(-(2q_yq_z - 2q_xq_w)) \\
\theta_y &=& \left\{ \begin{array}{ll}
  \arctan(\frac{2q_xq_z + 2q_yq_w}{2q_w^2 + 2q_z^2 - 1}) & (cos\theta_x \neq 0) \\
  \arctan(-\frac{2q_xq_z - 2q_yq_w}{2q_w^2 + 2q_x^2 - 1}) & (otherwise)
\end{array} \right. \\
\theta_z &=& \left\{ \begin{array}{ll}
  \arctan(\frac{2q_xq_y + 2q_zq_w}{2q_w^2 + 2q_y^2 - 1}) & (cos\theta_x \neq 0) \\
  0 & (otherwise)
\end{array} \right. \\
\end{eqnarray}
```

### 回転順YZX

```math
\begin{eqnarray}
\theta_x &=& \left\{ \begin{array}{ll}
  \arctan(-\frac{m_{12}}{m_{11}}) & (cos\theta_z \neq 0) \\
  0 & (otherwise)
\end{array} \right. \\
\theta_y &=& \left\{ \begin{array}{ll}
  \arctan(-\frac{m_{20}}{m_{00}}) & (cos\theta_z \neq 0) \\
  \arctan(\frac{m_{02}}{m_{22}}) & (otherwise)
\end{array} \right. \\
\theta_z &=& \arcsin(m_{10}) \\
\end{eqnarray}
```

より、

```math
\begin{eqnarray}
\theta_x &=& \left\{ \begin{array}{ll}
  \arctan(-\frac{2q_yq_z - 2q_xq_w}{2q_w^2 + 2q_y^2 - 1}) & (cos\theta_z \neq 0) \\
  0 & (otherwise)
\end{array} \right. \\
\theta_y &=& \left\{ \begin{array}{ll}
  \arctan(-\frac{2q_xq_z - 2q_yq_w}{2q_w^2 + 2q_x^2 - 1}) & (cos\theta_z \neq 0) \\
  \arctan(\frac{2q_xq_z + 2q_yq_w}{2q_w^2 + 2q_z^2 - 1}) & (otherwise)
\end{array} \right. \\
\theta_z &=& \arcsin(2q_xq_y + 2q_zq_w) \\
\end{eqnarray}
```

### 回転順ZXY

```math
\begin{eqnarray}
\theta_x &=& \arcsin(m_{21}) \\
\theta_y &=& \left\{ \begin{array}{ll}
  \arctan(-\frac{m_{20}}{m_{22}}) & (cos\theta_x \neq 0) \\
  0 & (otherwise)
\end{array} \right. \\
\theta_z &=& \left\{ \begin{array}{ll}
  \arctan(-\frac{m_{01}}{m_{11}}) & (cos\theta_x \neq 0) \\
  \arctan(\frac{m_{10}}{m_{00}}) & (otherwise)
\end{array} \right. \\
\end{eqnarray}
```

より、

```math
\begin{eqnarray}
\theta_x &=& \arcsin(2q_yq_z + 2q_xq_w) \\
\theta_y &=& \left\{ \begin{array}{ll}
  \arctan(-\frac{2q_xq_z - 2q_yq_w}{2q_w^2 + 2q_z^2 - 1}) & (cos\theta_x \neq 0) \\
  0 & (otherwise)
\end{array} \right. \\
\theta_z &=& \left\{ \begin{array}{ll}
  \arctan(-\frac{2q_xq_y - 2q_zq_w}{2q_w^2 + 2q_y^2 - 1}) & (cos\theta_x \neq 0) \\
  \arctan(\frac{2q_xq_y + 2q_zq_w}{2q_w^2 + 2q_x^2 - 1}) & (otherwise)
\end{array} \right. \\
\end{eqnarray}
```

### 回転順ZYX

```math
\begin{eqnarray}
\theta_x &=& \left\{ \begin{array}{ll}
  \arctan(\frac{m_{21}}{m_{22}}) & (cos\theta_y \neq 0) \\
  0 & (otherwise)
\end{array} \right. \\
\theta_y &=& \arcsin(-m_{20}) \\
\theta_z &=& \left\{ \begin{array}{ll}
  \arctan(\frac{m_{10}}{m_{00}}) & (cos\theta_y \neq 0) \\
  \arctan(-\frac{m_{01}}{m_{11}}) & (otherwise)
\end{array} \right. \\
\end{eqnarray}
```

より、

```math
\begin{eqnarray}
\theta_x &=& \left\{ \begin{array}{ll}
  \arctan(\frac{2q_yq_z + 2q_xq_w}{2q_w^2 + 2q_z^2 - 1}) & (cos\theta_y \neq 0) \\
  0 & (otherwise)
\end{array} \right. \\
\theta_y &=& \arcsin(-(2q_xq_z - 2q_yq_w)) \\
\theta_z &=& \left\{ \begin{array}{ll}
  \arctan(\frac{2q_xq_y + 2q_zq_w}{2q_w^2 + 2q_x^2 - 1}) & (cos\theta_y \neq 0) \\
  \arctan(-\frac{2q_xq_y - 2q_zq_w}{2q_w^2 + 2q_y^2 - 1}) & (otherwise)
\end{array} \right. \\
\end{eqnarray}
```

---

**追記**

実装編も書きました。
[回転行列、クォータニオン\(四元数\)、オイラー角の相互変換 \(実装編\) \- Qiita](https://qiita.com/aa_debdeb/items/abe90a9bd0b4809813da)
