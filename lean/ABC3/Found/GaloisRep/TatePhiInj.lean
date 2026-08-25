import ABC3.Found.GaloisRep.TatePhi

/-!
# Galois (G6) 第 285 ブロック —— **★★★★★★★★★★`Φ` は単射**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★到達点

> `Φ : Kˣ/q^ℤ → E_q(K)` は**単射**(`tatePhi_injective`)

★★★葉 (d)(単射性)が、3 領域すべてを含む形で閉じた。

## ★★★★★★★★★5 通りの貼り合わせ

`tate_inj_all` が芯である。局所環の二分法 `∀ x, IsUnit x ∨ x ∈ I` のもとで
`1 − a`・`1 − a'` と `a`・`a'` について場合分けする:

| `1 − a` | `1 − a'` | 根拠 |
|---|---|---|
| 非単元(原点近傍) | 非単元 | `tateZ_inj`(第 276) |
| 非単元 | 単元 | `tateXK_ne_of_origin`(第 280)——矛盾 |
| 単元 | 非単元 | 同上 |
| 単元、`a` 単元 | 単元、`a'` 単元 | `tate_point_inj_unit`(第 267) |
| 単元、`a ∈ I` | 単元、`a' ∈ I` | `tate_inj`(第 259) |
| 単元、片方だけ単元 | | `ne_of_isUnit_of_mem`(第 272)——矛盾 |

★★★★**「領域が違えば矛盾、同じなら領域内の単射性」**——これだけである。
矛盾の 3 通りは座標の住む場所(単元か、`I` の元か、`R` の外か)で決まる。

## ★★★局所環の二分法は仮定として受ける

`TateSetup` には入れず、`hloc : ∀ x : R, IsUnit x ∨ x ∈ I` として渡す。
★実際の適用では `R` は完備離散付値環、`I = 𝔪` なのでこれは成り立つ。

## ★★類の水準へ

`a(c) = a(c')` から `normRep c = normRep c'`、`mk ∘ normRep = id` で `c = c'`。
★単位類の場合は「原点に行くのは単位類だけ」(第 284)で処理する。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tate_inj_all` | ★★★★★★★★★**母数座標の単射性(全領域)** |
| `tatePhi_injective` | ★★★★★★★★★★**`Φ` は単射** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real WeierstrassCurve WeierstrassCurve.Affine QuotientGroup

/-! ## ★★★★★★★★★5 通りの貼り合わせ -/

section All

variable {R K : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  [Field K] [Algebra R K]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★**母数座標の単射性(全領域)**——5 通りの貼り合わせ。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_inj_all (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hinj : Function.Injective (algebraMap R K))
    (a w a' w' q : R) (hq : q ∈ I) (haw : a * w = q) (haw' : a' * w' = q)
    (hw : w ∈ I) (hw' : w' ∈ I)
    (hne : algebraMap R K (1 - a) ≠ 0) (hne' : algebraMap R K (1 - a') ≠ 0)
    (hX : tateXK a w q hq = (tateXK a' w' q hq : K))
    (hY : tateYK a w q hq = (tateYK a' w' q hq : K)) : a = a' := by
  rcases hloc (1 - a) with hu | hm
  · rcases hloc (1 - a') with hu' | hm'
    · have hXR : tateXpair a w q hq = tateXpair a' w' q hq := by
        refine hinj ?_
        rw [← tateXK_eq a w q hq hu, ← tateXK_eq a' w' q hq hu']
        exact hX
      have hYR : tateYpair a w q hq = tateYpair a' w' q hq := by
        refine hinj ?_
        rw [← tateYK_eq a w q hq hu, ← tateYK_eq a' w' q hq hu']
        exact hY
      rcases hloc a with hau | ham
      · rcases hloc a' with hau' | ham'
        · have e1 : wOf q a = w := wOf_eq_of_mul hau haw
          have e2 : wOf q a' = w' := wOf_eq_of_mul hau' haw'
          exact tate_point_inj_unit a a' q hq hau hau' hu hu'
            (by rw [e1, e2]; exact hXR) (by rw [e1, e2]; exact hYR)
        · exact absurd hXR (ne_of_isUnit_of_mem (isUnit_tateXpair a w q hq hw hau hu)
            (tateXpair_mem a' w' q hq ham' hw'))
      · rcases hloc a' with hau' | ham'
        · exact absurd hXR.symm (ne_of_isUnit_of_mem
            (isUnit_tateXpair a' w' q hq hw' hau' hu')
            (tateXpair_mem a w q hq ham hw))
        · exact (tate_inj a w a' w' q hq ham hw ham' hw' haw hXR hYR).1
    · exact absurd hX.symm (tateXK_ne_of_origin hinj a' w' a w q hq hm' hne' hu)
  · rcases hloc (1 - a') with hu' | hm'
    · exact absurd hX (tateXK_ne_of_origin hinj a w a' w' q hq hm hne hu')
    · have hau : IsUnit a := isUnit_of_one_sub_mem hm
      have hau' : IsUnit a' := isUnit_of_one_sub_mem hm'
      have e1 : wOf q a = w := wOf_eq_of_mul hau haw
      have e2 : wOf q a' = w' := wOf_eq_of_mul hau' haw'
      refine tateZ_inj a a' q hq hm hm' (hinj ?_)
      rw [tateZ_map a (wOf q a) q hq hm hne, tateZ_map a' (wOf q a') q hq hm' hne',
        e1, e2, hX, hY]

end All

/-! ## ★★★★★★★★★★`Φ` は単射 -/

section Phi

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [DecidableEq K] [Algebra R K]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**`Φ` は単射**——葉 (d) が 3 領域すべてを含む形で閉じた。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tatePhi_injective (S : TateSetup R I K) (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0) :
    Function.Injective (tatePhi S hΔ) := by
  intro c c' h
  by_cases hc : c = 1
  · subst hc
    rw [tatePhi_one] at h
    by_contra hne
    exact tatePhi_ne_zero S hΔ (fun hh => hne hh.symm) h.symm
  · by_cases hc' : c' = 1
    · subst hc'
      rw [tatePhi_one] at h
      exact absurd h (tatePhi_ne_zero S hΔ hc)
    · rw [tatePhi_eq S hΔ hc, tatePhi_eq S hΔ hc', tatePtPair, tatePtPair] at h
      simp only [Point.some.injEq] at h
      have haa : tateAOf S c = tateAOf S c' :=
        tate_inj_all hloc S.hinj _ _ _ _ S.q S.hq (tateAOf_mul_tateWOf S c)
          (tateAOf_mul_tateWOf S c') (tateWOf_mem S c) (tateWOf_mem S c')
          (tateAOf_ne_one S hc) (tateAOf_ne_one S hc') h.1 h.2
      have hnr : (normRep S.v S.Q S.hQ c : K) = (normRep S.v S.Q S.hQ c' : K) := by
        rw [← (tateAOf_spec S c).1, ← (tateAOf_spec S c').1, haa]
      have h2 : normRep S.v S.Q S.hQ c = normRep S.v S.Q S.hQ c' := Units.ext hnr
      rw [← normRep_mk S.v S.Q S.hQ c, ← normRep_mk S.v S.Q S.hQ c', h2]

end Phi

/-! ## ★出典の紐付け(`.src`) -/

def tate_inj_all.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——母数座標の単射性(全領域))",
    sectionId := "genell-def-3-3" }

def tatePhi_injective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——Phi は単射)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
