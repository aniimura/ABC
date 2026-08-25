import ABC3.Found.GaloisRep.TateSurjUnit

/-!
# Galois (G6) 第 267 ブロック —— **★★★★★★★★単元の領域での単射性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★母数座標は単射

第 259 の単射性は `a, w ∈ I`(**環帯**の領域)での逐次近似だった。そこでは
`X ≡ a + w`、`Y ≡ −w` という主要項を使ったが、★`a` が**単元**のときはこの展開が
効かない(`f(a) − a` が位を上げない)。

★★★しかし第 265 の `κ = Λ − id` の評価がそのまま使える:

    Λ(a) = Λ(b)  かつ  a − b ∈ I^j   ⟹   a − b = κ(b) − κ(a) ∈ I^{j+1}

を繰り返して `IsHausdorff` で `a = b`(`tate_lambda_inj`)。
★★**同じ縮小性が全射性(第 266)と単射性(本ブロック)の両方を出す**——
不動点の存在と一意性が 1 つの評価から出るのは、縮小写像定理の常である。

## ★★`X` だけの版

曲線の式の差から `(Y − Y')(Y + Y' + X) = 0`。`R` が整域なので分岐は 2 つで、
第 2 の枝は `P(b) = −P(a)` を意味する(`tate_inj_X_unit`)。

## ★★領域ごとの単射性がそろった

| 領域 | 定理 |
|---|---|
| 環帯(`a, w ∈ I`) | `tate_inj`・`tate_inj_X`(第 259) |
| 単元(`a` 単元、`w ∈ I`) | `tate_lambda_inj`・`tate_inj_X_unit`(本ブロック) |

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tate_lambda_inj` | ★★★★★★★★**母数座標は単射** |
| `tate_point_inj_unit` | ★★★★★★★★座標の単射性(単元の領域) |
| `tate_inj_X_unit` | ★★★★★★★★`X` だけの版 |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★★★★★★★★**母数座標は単射**——単元の領域での葉 (d)。

★第 265 の `κ` の評価をそのまま逐次近似に使う。全射性(第 266)と同じ縮小性から出る。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_lambda_inj [IsAdicComplete I R] (a b q : R) (hq : q ∈ I)
    (hau : IsUnit a) (hbu : IsUnit b) (ha1 : IsUnit (1 - a)) (hb1 : IsUnit (1 - b))
    (hΛ : tateLambda a (wOf q a) q hq = tateLambda b (wOf q b) q hq) : a = b := by
  have hstep : ∀ j : ℕ, a - b ∈ I ^ j := by
    intro j
    induction j with
    | zero => simp
    | succ n ih =>
      have h := tateKappa_diff_mem a b q hq hau hbu ha1 hb1 ih
      have heq : (tateLambda a (wOf q a) q hq - a) - (tateLambda b (wOf q b) q hq - b)
          = -(a - b) := by rw [hΛ]; ring
      rw [heq] at h
      have h2 := neg_mem h
      rwa [neg_neg] at h2
  exact sub_eq_zero.1 (eq_zero_of_mem_pow (I := I) hstep)

/-- ★★★★★★★★**単元の領域での座標の単射性**。 -/
theorem tate_point_inj_unit [IsAdicComplete I R] (a b q : R) (hq : q ∈ I)
    (hau : IsUnit a) (hbu : IsUnit b) (ha1 : IsUnit (1 - a)) (hb1 : IsUnit (1 - b))
    (hX : tateXpair a (wOf q a) q hq = tateXpair b (wOf q b) q hq)
    (hY : tateYpair a (wOf q a) q hq = tateYpair b (wOf q b) q hq) : a = b := by
  refine tate_lambda_inj a b q hq hau hbu ha1 hb1 ?_
  rw [tateLambda, tateLambda, hX, hY]

/-- ★★★★★★★★**`X` だけでも `±` を除いて決まる**(単元の領域)。 -/
theorem tate_inj_X_unit [IsAdicComplete I R] [IsDomain R] (a b q : R) (hq : q ∈ I)
    (hau : IsUnit a) (hbu : IsUnit b) (ha1 : IsUnit (1 - a)) (hb1 : IsUnit (1 - b))
    (hX : tateXpair a (wOf q a) q hq = tateXpair b (wOf q b) q hq) :
    a = b ∨ tateYpair b (wOf q b) q hq
      = -tateXpair a (wOf q a) q hq - tateYpair a (wOf q a) q hq := by
  have e1 := tate_equation a (wOf q a) q hq (mul_wOf q a hau) ha1
    (isUnit_one_sub (wOf_mem a hq))
  have e2 := tate_equation b (wOf q b) q hq (mul_wOf q b hbu) hb1
    (isUnit_one_sub (wOf_mem b hq))
  rw [← hX] at e2
  have hfac : (tateYpair a (wOf q a) q hq - tateYpair b (wOf q b) q hq)
      * (tateYpair a (wOf q a) q hq + tateYpair b (wOf q b) q hq
        + tateXpair a (wOf q a) q hq) = 0 := by
    linear_combination e1 - e2
  rcases mul_eq_zero.1 hfac with h | h
  · exact Or.inl (tate_point_inj_unit a b q hq hau hbu ha1 hb1 hX (sub_eq_zero.1 h))
  · exact Or.inr (by linear_combination h)

/-! ## ★出典の紐付け(`.src`) -/

def tate_lambda_inj.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——母数座標は単射)",
    sectionId := "genell-def-3-3" }

def tate_inj_X_unit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——単元の領域で X だけで ± を除いて決まる)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
