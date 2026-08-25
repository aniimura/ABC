import ABC3.Found.GaloisRep.TateFormal

/-!
# Galois (G6) 第 277 ブロック —— **★★★★★★分母を払った共線性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★なぜ共線性にも分母払いが要るか

葉 (e) の 3 領域はそろった(第 276)が、**群法則はまだ原点近傍を含まない**。
`tate_collinear`(第 256)は `hcp : ∀ i, IsUnit (1 − collPts u v w i)` を要求し、
`u ≡ 1` ではこれが破れるからである。

★★★退化した場合(倍化)を補助母数で処理するには、補助母数を**どの領域からでも**
取れねばならない。剰余体が `𝔽₂` で `v(q) = 1` のときは、単元でも環帯でもない元しか
無い——**原点近傍を含まない群法則は空虚になりうる**。

## ★★★★★★行ごとに分母を払う

`collDefect` は行列式

    det [[X₁, Y₁, 1], [X₂, Y₂, 1], [X₃, Y₃, 1]]

である。★★★**行 `i` を `M_i` 倍すれば行列式は `M₁M₂M₃` 倍になる**から、

    M_i := (1−p_i)³(1−r_i)³        (`p_i` は母数、`r_i` は相方)

と置いて `XE_i := M_i·X_i`、`YE_i := M_i·Y_i`、そして第 3 列を `M_i` にすればよい。
★★これで **`Ring.inverse` を一切含まない行列式**になる。

## ★★★★★★両側の分母を払う

第 274 は `1 − a` の側だけを払ったが、共線性では相方 `r` の側も要る:

    XE = (1−p)(1−r)³·p + (1−p)³(1−r)·r + M·(尾)
    YE = (1−r)³·p² − (1−p)³(1−r)·r − (1−p)³·r² + M·(尾)

★`(1−p)²f(p) = p`、`(1−p)³g(p) = p²`(第 262)を**両方の変数に**当てるだけである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tateM` | ★★★対の分母 `(1−p)³(1−r)³` |
| `tateXpairEE`・`tateYpairEE` | ★★★★★★両側の分母を払った座標 |
| `tateXtruncEE`・`tateYtruncEE` | ★★★★★★その切り詰め |
| `collDefectE`・`collDefectTruncE` | ★★★★★★**分母を払った共線性の差** |
| `collDefectE_eq` | ★★★★★`collDefectE = M₁M₂M₃·collDefect` |
| `collDefectE_sub_trunc` | ★★★★差と切り詰めの差は `I^n` |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★★★★両側の分母を払った座標 -/

/-- ★★★**対の分母**——`(1−p)³(1−r)³`。 -/
noncomputable def tateM (p r : R) : R := (1 - p) ^ 3 * (1 - r) ^ 3

/-- ★★★★★★**両側の分母を払った `X`**。 -/
noncomputable def tateXpairEE [IsAdicComplete I R] (p r q : R) (hq : q ∈ I) : R :=
  (1 - p) * (1 - r) ^ 3 * p + (1 - p) ^ 3 * (1 - r) * r
    + tateM p r * (tateXtail p q hq + tateXtail r q hq - 2 * evalAdic (sigmaSeries 1) q hq)

/-- ★★★★★★**両側の分母を払った `Y`**。 -/
noncomputable def tateYpairEE [IsAdicComplete I R] (p r q : R) (hq : q ∈ I) : R :=
  (1 - r) ^ 3 * p ^ 2 - (1 - p) ^ 3 * (1 - r) * r - (1 - p) ^ 3 * r ^ 2
    + tateM p r * (tateYtail p q hq - tateXtail r q hq - tateYtail r q hq
      + evalAdic (sigmaSeries 1) q hq)

theorem tateXpairEE_eq [IsAdicComplete I R] (p r q : R) (hq : q ∈ I)
    (hp : IsUnit (1 - p)) (hr : IsUnit (1 - r)) :
    tateXpairEE p r q hq = tateM p r * tateXpair p r q hq := by
  rw [tateXpairEE, tateXpair, tateM]
  linear_combination (-((1 - p) * (1 - r) ^ 3)) * mul_tateXterm' hp
    + (-((1 - p) ^ 3 * (1 - r))) * mul_tateXterm' hr

theorem tateYpairEE_eq [IsAdicComplete I R] (p r q : R) (hq : q ∈ I)
    (hp : IsUnit (1 - p)) (hr : IsUnit (1 - r)) :
    tateYpairEE p r q hq = tateM p r * tateYpair p r q hq := by
  rw [tateYpairEE, tateYpair, tateM]
  linear_combination (-((1 - r) ^ 3)) * mul_tateYterm' hp
    + ((1 - p) ^ 3 * (1 - r)) * mul_tateXterm' hr
    + ((1 - p) ^ 3) * mul_tateYterm' hr

/-! ## ★★★★★★その切り詰め -/

noncomputable def tateXtruncEE (n : ℕ) (p r q : R) : R :=
  (1 - p) * (1 - r) ^ 3 * p + (1 - p) ^ 3 * (1 - r) * r
    + tateM p r * (partialSum (fun m => tateXterm (q ^ (m + 1) * p)) n
      + partialSum (fun m => tateXterm (q ^ (m + 1) * r)) n
      - 2 * partialEval (sigmaSeries 1) q n)

noncomputable def tateYtruncEE (n : ℕ) (p r q : R) : R :=
  (1 - r) ^ 3 * p ^ 2 - (1 - p) ^ 3 * (1 - r) * r - (1 - p) ^ 3 * r ^ 2
    + tateM p r * (partialSum (fun m => tateYterm (q ^ (m + 1) * p)) n
      - partialSum (fun m => tateXterm (q ^ (m + 1) * r)) n
      - partialSum (fun m => tateYterm (q ^ (m + 1) * r)) n
      + partialEval (sigmaSeries 1) q n)

theorem tateXtruncEE_eq (n : ℕ) (p r q : R) (hp : IsUnit (1 - p)) (hr : IsUnit (1 - r)) :
    tateXtruncEE n p r q = tateM p r * tateXtrunc n p r q := by
  rw [tateXtruncEE, tateXtrunc, tateM]
  linear_combination (-((1 - p) * (1 - r) ^ 3)) * mul_tateXterm' hp
    + (-((1 - p) ^ 3 * (1 - r))) * mul_tateXterm' hr

theorem tateYtruncEE_eq (n : ℕ) (p r q : R) (hp : IsUnit (1 - p)) (hr : IsUnit (1 - r)) :
    tateYtruncEE n p r q = tateM p r * tateYtrunc n p r q := by
  rw [tateYtruncEE, tateYtrunc, tateM]
  linear_combination (-((1 - r) ^ 3)) * mul_tateYterm' hp
    + ((1 - p) ^ 3 * (1 - r)) * mul_tateXterm' hr
    + ((1 - p) ^ 3) * mul_tateYterm' hr

theorem tateXpairEE_sub_trunc [IsAdicComplete I R] (n : ℕ) (p r q : R) (hq : q ∈ I) :
    tateXpairEE p r q hq - tateXtruncEE n p r q ∈ I ^ n := by
  have hkey : tateXpairEE p r q hq - tateXtruncEE n p r q
      = tateM p r * ((tateXtail p q hq - partialSum (fun m => tateXterm (q ^ (m + 1) * p)) n)
        + (tateXtail r q hq - partialSum (fun m => tateXterm (q ^ (m + 1) * r)) n)
        - 2 * (evalAdic (sigmaSeries 1) q hq - partialEval (sigmaSeries 1) q n)) := by
    rw [tateXpairEE, tateXtruncEE]; ring
  rw [hkey]
  exact Ideal.mul_mem_left _ _ (Ideal.sub_mem _ (Ideal.add_mem _
    (tateXtail_sub_partialSum n p q hq) (tateXtail_sub_partialSum n r q hq))
    (Ideal.mul_mem_left _ _ (evalAdic_sub_partialEval n (sigmaSeries 1) q hq)))

theorem tateYpairEE_sub_trunc [IsAdicComplete I R] (n : ℕ) (p r q : R) (hq : q ∈ I) :
    tateYpairEE p r q hq - tateYtruncEE n p r q ∈ I ^ n := by
  have hkey : tateYpairEE p r q hq - tateYtruncEE n p r q
      = tateM p r * ((tateYtail p q hq - partialSum (fun m => tateYterm (q ^ (m + 1) * p)) n)
        - (tateXtail r q hq - partialSum (fun m => tateXterm (q ^ (m + 1) * r)) n)
        - (tateYtail r q hq - partialSum (fun m => tateYterm (q ^ (m + 1) * r)) n)
        + (evalAdic (sigmaSeries 1) q hq - partialEval (sigmaSeries 1) q n)) := by
    rw [tateYpairEE, tateYtruncEE]; ring
  rw [hkey]
  exact Ideal.mul_mem_left _ _ (Ideal.add_mem _ (Ideal.sub_mem _ (Ideal.sub_mem _
    (tateYtail_sub_partialSum n p q hq) (tateXtail_sub_partialSum n r q hq))
    (tateYtail_sub_partialSum n r q hq)) (evalAdic_sub_partialEval n (sigmaSeries 1) q hq))

/-! ## ★★★★★★分母を払った共線性の差 -/

/-- ★★★★★★**分母を払った共線性の差**——行 `i` を `M_i` 倍した行列式。 -/
noncomputable def collDefectE [IsAdicComplete I R] (u v w q : R) (hq : q ∈ I) : R :=
  tateXpairEE u (v * w) q hq * (tateYpairEE v (u * w) q hq * tateM w (u * v)
      - tateYpairEE w (u * v) q hq * tateM v (u * w))
    + tateXpairEE v (u * w) q hq * (tateYpairEE w (u * v) q hq * tateM u (v * w)
      - tateYpairEE u (v * w) q hq * tateM w (u * v))
    + tateXpairEE w (u * v) q hq * (tateYpairEE u (v * w) q hq * tateM v (u * w)
      - tateYpairEE v (u * w) q hq * tateM u (v * w))

/-- ★★★★★★分母を払った共線性の差の切り詰め。 -/
noncomputable def collDefectTruncE (n : ℕ) (u v w q : R) : R :=
  tateXtruncEE n u (v * w) q * (tateYtruncEE n v (u * w) q * tateM w (u * v)
      - tateYtruncEE n w (u * v) q * tateM v (u * w))
    + tateXtruncEE n v (u * w) q * (tateYtruncEE n w (u * v) q * tateM u (v * w)
      - tateYtruncEE n u (v * w) q * tateM w (u * v))
    + tateXtruncEE n w (u * v) q * (tateYtruncEE n u (v * w) q * tateM v (u * w)
      - tateYtruncEE n v (u * w) q * tateM u (v * w))

set_option maxHeartbeats 1600000 in
/-- ★★★★★**`collDefectE = M₁M₂M₃·collDefect`**——行を倍したぶんだけ行列式が倍になる。 -/
theorem collDefectE_eq [IsAdicComplete I R] (u v w q : R) (hq : q ∈ I)
    (hcp : ∀ i, IsUnit (1 - collPts u v w i)) :
    collDefectE u v w q hq
      = tateM u (v * w) * tateM v (u * w) * tateM w (u * v) * collDefect u v w q hq := by
  have h0 : IsUnit (1 - u) := hcp 0
  have h1 : IsUnit (1 - v) := hcp 1
  have h2 : IsUnit (1 - w) := hcp 2
  have h3 : IsUnit (1 - v * w) := hcp 3
  have h4 : IsUnit (1 - u * w) := hcp 4
  have h5 : IsUnit (1 - u * v) := hcp 5
  rw [collDefectE, collDefect,
    tateXpairEE_eq u (v * w) q hq h0 h3, tateYpairEE_eq u (v * w) q hq h0 h3,
    tateXpairEE_eq v (u * w) q hq h1 h4, tateYpairEE_eq v (u * w) q hq h1 h4,
    tateXpairEE_eq w (u * v) q hq h2 h5, tateYpairEE_eq w (u * v) q hq h2 h5]
  ring

set_option maxHeartbeats 1600000 in
theorem collDefectTruncE_eq (n : ℕ) (u v w q : R)
    (hcp : ∀ i, IsUnit (1 - collPts u v w i)) :
    collDefectTruncE n u v w q
      = tateM u (v * w) * tateM v (u * w) * tateM w (u * v) * collDefectTrunc n u v w q := by
  have h0 : IsUnit (1 - u) := hcp 0
  have h1 : IsUnit (1 - v) := hcp 1
  have h2 : IsUnit (1 - w) := hcp 2
  have h3 : IsUnit (1 - v * w) := hcp 3
  have h4 : IsUnit (1 - u * w) := hcp 4
  have h5 : IsUnit (1 - u * v) := hcp 5
  rw [collDefectTruncE, collDefectTrunc,
    tateXtruncEE_eq n u (v * w) q h0 h3, tateYtruncEE_eq n u (v * w) q h0 h3,
    tateXtruncEE_eq n v (u * w) q h1 h4, tateYtruncEE_eq n v (u * w) q h1 h4,
    tateXtruncEE_eq n w (u * v) q h2 h5, tateYtruncEE_eq n w (u * v) q h2 h5]
  ring

set_option maxHeartbeats 1600000 in
theorem collDefectE_sub_trunc [IsAdicComplete I R] (n : ℕ) (u v w q : R) (hq : q ∈ I) :
    collDefectE u v w q hq - collDefectTruncE n u v w q ∈ I ^ n := by
  refine Ideal.Quotient.eq.1 ?_
  have hX1 : (Ideal.Quotient.mk (I ^ n)) (tateXpairEE u (v * w) q hq)
      = (Ideal.Quotient.mk (I ^ n)) (tateXtruncEE n u (v * w) q) :=
    Ideal.Quotient.eq.2 (tateXpairEE_sub_trunc n u (v * w) q hq)
  have hX2 : (Ideal.Quotient.mk (I ^ n)) (tateXpairEE v (u * w) q hq)
      = (Ideal.Quotient.mk (I ^ n)) (tateXtruncEE n v (u * w) q) :=
    Ideal.Quotient.eq.2 (tateXpairEE_sub_trunc n v (u * w) q hq)
  have hX3 : (Ideal.Quotient.mk (I ^ n)) (tateXpairEE w (u * v) q hq)
      = (Ideal.Quotient.mk (I ^ n)) (tateXtruncEE n w (u * v) q) :=
    Ideal.Quotient.eq.2 (tateXpairEE_sub_trunc n w (u * v) q hq)
  have hY1 : (Ideal.Quotient.mk (I ^ n)) (tateYpairEE u (v * w) q hq)
      = (Ideal.Quotient.mk (I ^ n)) (tateYtruncEE n u (v * w) q) :=
    Ideal.Quotient.eq.2 (tateYpairEE_sub_trunc n u (v * w) q hq)
  have hY2 : (Ideal.Quotient.mk (I ^ n)) (tateYpairEE v (u * w) q hq)
      = (Ideal.Quotient.mk (I ^ n)) (tateYtruncEE n v (u * w) q) :=
    Ideal.Quotient.eq.2 (tateYpairEE_sub_trunc n v (u * w) q hq)
  have hY3 : (Ideal.Quotient.mk (I ^ n)) (tateYpairEE w (u * v) q hq)
      = (Ideal.Quotient.mk (I ^ n)) (tateYtruncEE n w (u * v) q) :=
    Ideal.Quotient.eq.2 (tateYpairEE_sub_trunc n w (u * v) q hq)
  simp only [collDefectE, collDefectTruncE, map_sub, map_add, map_mul, hX1, hX2, hX3,
    hY1, hY2, hY3]

/-! ## ★出典の紐付け(`.src`) -/

def collDefectE.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——分母を払った共線性の差)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
