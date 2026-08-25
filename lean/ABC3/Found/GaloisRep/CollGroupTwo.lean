import ABC3.Found.GaloisRep.TatePhiTwo

/-!
# Galois (G6) 第 295 ブロック —— **★★★★★★★★★★単元 2 つの場合の共線性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★なぜ `n = 0` が最後まで残るのか

`v(q) = 1` のときは**すべての類の付値が 0**なので、`c₁c₂c₃ = 1` なら必ず `n = 0` である。
★★★つまり `n = 0` は端の場合ではなく、**`v(q) = 1` では唯一の場合**である。
逆元で `n = 1` に落とす手(第 294)も使えない。正面から扱うほかない。

## ★★★★★★★★向きを入れ替えれば座標が定義できる

`n = 0` では `a₁a₂a₃ = 1`(3 つとも単元)なので、3 つ組 `(a₁, a₂, q·a₃)` を取ると
相方は `q/a₁`・`q/a₂`・`a₁a₂` である。★最初の 2 つは `𝔪` に入るが、
**`a₁a₂ = a₃⁻¹` は単元**なので、第 3 の点の座標 `X(q a₃, a₁a₂)` は
`1 − a₁a₂` が単元でないと定義できない(`a₃ ≡ 1` のとき破れる)。

★★★★**向きを入れ替える**:第 3 の点を `(a₁a₂, q a₃)` の順で書けば、母数が `a₁a₂`、
相方が `q a₃ ∈ 𝔪` になり、**`1 − 相方` は常に単元**である。

## ★★★★★両側払いは対称、`Y` は反転する

    XE(r,p) = XE(p,r)                            (`tateXpairEE_symm`)
    YE(r,p) = −XE(p,r) − YE(p,r)                 (`tateYpairEE_swap`)

★★どちらも**仮定なしの多項式恒等式**である(分母を払ってあるから `ring` で出る)。
★★★これで第 278 の `collDefectE_eq_zero`(仮定なし)の第 3 行を、
入れ替えた向きの座標で書き直せる。

## ★★★★★★★★★★到達点

> `u·v·w = q`、`vw, uw, w ∈ 𝔪`、`1 − u, 1 − v, 1 − uv` が `K` で `0` でない —— このとき
> `X₁(Y₂ − Y₃) + X₂(Y₃ − Y₁) + X₃(Y₁ − Y₂) = 0`(`Y₃ := −X₃'' − Y₃''`)

★`(X₃'', Y₃'')` は入れ替えた向きの座標、`Y₃` はもとの向きの `Y` にあたる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tateXpairEE_symm`・`tateYpairEE_swap`・`tateM_symm` | ★★★★★両側払いの対称性 |
| `collDet_K_eq_zero_two` | ★★★★★★★★★★**単元 2 つの場合の共線性** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real WeierstrassCurve WeierstrassCurve.Affine

/-! ## ★★★★★両側払いの対称性 -/

section Symm

variable {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]

/-- ★★★★★**両側払いの `X` は対称**。 -/
theorem tateXpairEE_symm (p r q : R) (hq : q ∈ I) :
    tateXpairEE r p q hq = tateXpairEE p r q hq := by
  rw [tateXpairEE, tateXpairEE, tateM, tateM]
  ring

/-- ★★★★★**両側払いの `Y` の交換則**——`YE(r,p) = −XE(p,r) − YE(p,r)`。 -/
theorem tateYpairEE_swap (p r q : R) (hq : q ∈ I) :
    tateYpairEE r p q hq = -tateXpairEE p r q hq - tateYpairEE p r q hq := by
  rw [tateYpairEE, tateYpairEE, tateXpairEE, tateM, tateM]
  ring

theorem tateM_symm (p r : R) : tateM r p = tateM p r := by
  rw [tateM, tateM]; ring

end Symm

/-! ## ★★★★★★★★★★単元 2 つの場合の共線性 -/

section Two

variable {R K : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R] [Field K] [Algebra R K]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**`K` の水準の共線性(単元 2 つの場合)**——第 3 の点を入れ替えた向きで。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem collDet_K_eq_zero_two (u v w q : R) (hq : q ∈ I) (huvw : u * v * w = q)
    (hvw : v * w ∈ I) (huw : u * w ∈ I) (hw : w ∈ I)
    (hu : algebraMap R K (1 - u) ≠ 0) (hv : algebraMap R K (1 - v) ≠ 0)
    (huv : algebraMap R K (1 - u * v) ≠ 0) :
    tateXK u (v * w) q hq
        * (tateYK (K := K) v (u * w) q hq
          - (-tateXK (u * v) w q hq - tateYK (K := K) (u * v) w q hq))
      + tateXK v (u * w) q hq
        * ((-tateXK (u * v) w q hq - tateYK (K := K) (u * v) w q hq)
          - tateYK (K := K) u (v * w) q hq)
      + tateXK (u * v) w q hq
        * (tateYK (K := K) u (v * w) q hq - tateYK (K := K) v (u * w) q hq) = 0 := by
  have uvw : IsUnit (1 - v * w) := isUnit_one_sub hvw
  have uuw : IsUnit (1 - u * w) := isUnit_one_sub huw
  have uw' : IsUnit (1 - w) := isUnit_one_sub hw
  have nvw : algebraMap R K (1 - v * w) ≠ 0 := (uvw.map (algebraMap R K)).ne_zero
  have nuw : algebraMap R K (1 - u * w) ≠ 0 := (uuw.map (algebraMap R K)).ne_zero
  have nw : algebraMap R K (1 - w) ≠ 0 := (uw'.map (algebraMap R K)).ne_zero
  have h0 := congrArg (algebraMap R K) (collDefectE_eq_zero u v w q hq huvw)
  rw [collDefectE, map_zero] at h0
  rw [tateXpairEE_symm (u * v) w q hq, tateYpairEE_swap (u * v) w q hq,
    tateM_symm (u * v) w] at h0
  simp only [map_add, map_sub, map_mul, map_neg] at h0
  rw [map_tateXpairEE u (v * w) q hq uvw hu, map_tateYpairEE u (v * w) q hq uvw hu,
    map_tateXpairEE v (u * w) q hq uuw hv, map_tateYpairEE v (u * w) q hq uuw hv,
    map_tateXpairEE (u * v) w q hq uw' huv, map_tateYpairEE (u * v) w q hq uw' huv] at h0
  have hM : (algebraMap R K (tateM u (v * w)) * algebraMap R K (tateM v (u * w))
      * algebraMap R K (tateM (u * v) w)) ≠ 0 :=
    mul_ne_zero (mul_ne_zero (map_tateM_ne_zero u (v * w) hu nvw)
      (map_tateM_ne_zero v (u * w) hv nuw)) (map_tateM_ne_zero (u * v) w huv nw)
  have hkey : (algebraMap R K (tateM u (v * w)) * algebraMap R K (tateM v (u * w))
      * algebraMap R K (tateM (u * v) w))
      * (tateXK u (v * w) q hq
          * (tateYK (K := K) v (u * w) q hq
            - (-tateXK (u * v) w q hq - tateYK (K := K) (u * v) w q hq))
        + tateXK v (u * w) q hq
          * ((-tateXK (u * v) w q hq - tateYK (K := K) (u * v) w q hq)
            - tateYK (K := K) u (v * w) q hq)
        + tateXK (u * v) w q hq
          * (tateYK (K := K) u (v * w) q hq - tateYK (K := K) v (u * w) q hq)) = 0 := by
    linear_combination h0
  exact (mul_eq_zero.1 hkey).resolve_left hM

end Two

/-! ## ★出典の紐付け(`.src`) -/

def tateYpairEE_swap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——両側払いの Y の交換則)",
    sectionId := "genell-def-3-3" }

def collDet_K_eq_zero_two.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——単元 2 つの場合の共線性)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
