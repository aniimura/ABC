import ABC3.Found.GaloisRep.TatePt

/-!
# Galois (G6) 第 283 ブロック —— **★★★★★★★点は類だけで決まる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★到達点——`Φ` が定まる根拠

> 同じ類 `c ∈ Kˣ/q^ℤ` から来る 2 つの正規化した対は同じ点を与える
> (`tatePtPair_eq_of_same_class'`)

★★★第 213 の `pair_eq_of_same_class` は「対 `(a,w)` そのものが一致する」を言うので、
**対の任意の関数が類だけで決まる**。第 282 の `tatePtPair` はまさにその関数である。

## ★★★★対の一意性の骨

| 段 | 内容 |
|---|---|
| 1 | 同じ類の正規化元は代表元に一致する(`eq_normRep`、第 213) |
| 2 | `algebraMap a = u = u' = algebraMap a'` と単射性で `a = a'` |
| 3 | `a·w = q = a'·w'` と `a ≠ 0` の相殺で `w = w'` |

★`a ≠ 0` は `algebraMap a = (u : K)` が**単元**だから出る。

## ★★★単位類の代表元は `1`

`vAdd v 1 = 0` は `[0, vAdd v q)` に入り、類も単位元なので `normRep v Q hQ 1 = 1`。
★★★したがって **`a = 1` ⟺ 類が単位元**であり、`Φ` の場合分けはここで切れる
(`1 − a = 0` になるのは単位類のときだけ)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `pair_eq_of_same_class'` | ★★★★★対は類だけで決まる(`q` 固定) |
| `tateXK_eq_of_same_class`・`tateYK_eq_of_same_class` | ★★★★★★座標も類だけで決まる |
| `tatePtPair_eq_of_same_class'` | ★★★★★★★**点も類だけで決まる** |
| `vAdd_one`・`normRep_one` | ★★★単位類の代表元は `1` |
-/

namespace ABC3.Found.GaloisRep

open Complex Real WeierstrassCurve WeierstrassCurve.Affine QuotientGroup

/-! ## ★★★単位類の代表元は `1` -/

section Rep

variable {K : Type} [Field K]

theorem vAdd_one (v : Kˣ →* Multiplicative ℤ) : vAdd v 1 = 0 := by
  simp [vAdd]

/-- ★★★**単位類の正規化代表元は `1`**。 -/
theorem normRep_one (v : Kˣ →* Multiplicative ℤ) (Q : Kˣ) (hQ : 0 < vAdd v Q) :
    normRep v Q hQ 1 = 1 := by
  refine (eq_normRep v Q hQ 1 1 ?_ ?_ ?_).symm
  · simp
  · rw [vAdd_one]
  · rw [vAdd_one]; exact hQ

end Rep

/-! ## ★★★★★★★点は類だけで決まる -/

section Class

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [Algebra R K]

/-- ★★★★★**対は類だけで決まる**(`q` を固定した形)。 -/
theorem pair_eq_of_same_class' (hinj : Function.Injective (algebraMap R K))
    (v : Kˣ →* Multiplicative ℤ) (Q : Kˣ) (hQ : 0 < vAdd v Q) (q : R)
    (hQq : algebraMap R K q = (Q : K))
    {a w a' w' : R} {u u' : Kˣ}
    (hu : algebraMap R K a = (u : K)) (hu' : algebraMap R K a' = (u' : K))
    (haw : a * w = q) (haw' : a' * w' = q)
    (h0 : 0 ≤ vAdd v u) (h1 : vAdd v u < vAdd v Q)
    (h0' : 0 ≤ vAdd v u') (h1' : vAdd v u' < vAdd v Q)
    (hcls : QuotientGroup.mk (s := Subgroup.zpowers Q) u = QuotientGroup.mk u') :
    a = a' ∧ w = w' :=
  pair_eq_of_same_class hinj v Q hQ hu hu' (by rw [haw, hQq]) (by rw [haw', hQq])
    h0 h1 h0' h1' hcls

/-- ★★★★★★**`K` の水準の `X` も類だけで決まる**。 -/
theorem tateXK_eq_of_same_class (hinj : Function.Injective (algebraMap R K))
    (v : Kˣ →* Multiplicative ℤ) (Q : Kˣ) (hQ : 0 < vAdd v Q)
    {a w a' w' : R} (hw : w ∈ I) (hw' : w' ∈ I) {u u' : Kˣ}
    (hu : algebraMap R K a = (u : K)) (hu' : algebraMap R K a' = (u' : K))
    (hqa : algebraMap R K (a * w) = (Q : K)) (hqa' : algebraMap R K (a' * w') = (Q : K))
    (h0 : 0 ≤ vAdd v u) (h1 : vAdd v u < vAdd v Q)
    (h0' : 0 ≤ vAdd v u') (h1' : vAdd v u' < vAdd v Q)
    (hcls : QuotientGroup.mk (s := Subgroup.zpowers Q) u = QuotientGroup.mk u') :
    (tateXK (I := I) a w (a * w) (Ideal.mul_mem_left _ _ hw) : K)
      = tateXK (I := I) a' w' (a' * w') (Ideal.mul_mem_left _ _ hw') := by
  obtain ⟨rfl, rfl⟩ := pair_eq_of_same_class hinj v Q hQ hu hu' hqa hqa' h0 h1 h0' h1' hcls
  rfl

theorem tateYK_eq_of_same_class (hinj : Function.Injective (algebraMap R K))
    (v : Kˣ →* Multiplicative ℤ) (Q : Kˣ) (hQ : 0 < vAdd v Q)
    {a w a' w' : R} (hw : w ∈ I) (hw' : w' ∈ I) {u u' : Kˣ}
    (hu : algebraMap R K a = (u : K)) (hu' : algebraMap R K a' = (u' : K))
    (hqa : algebraMap R K (a * w) = (Q : K)) (hqa' : algebraMap R K (a' * w') = (Q : K))
    (h0 : 0 ≤ vAdd v u) (h1 : vAdd v u < vAdd v Q)
    (h0' : 0 ≤ vAdd v u') (h1' : vAdd v u' < vAdd v Q)
    (hcls : QuotientGroup.mk (s := Subgroup.zpowers Q) u = QuotientGroup.mk u') :
    (tateYK (I := I) a w (a * w) (Ideal.mul_mem_left _ _ hw) : K)
      = tateYK (I := I) a' w' (a' * w') (Ideal.mul_mem_left _ _ hw') := by
  obtain ⟨rfl, rfl⟩ := pair_eq_of_same_class hinj v Q hQ hu hu' hqa hqa' h0 h1 h0' h1' hcls
  rfl

/-- ★★★★★★★**点も類だけで決まる**——`Φ` が定まるための核心。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tatePtPair_eq_of_same_class' [DecidableEq K]
    (hinj : Function.Injective (algebraMap R K))
    (v : Kˣ →* Multiplicative ℤ) (Q : Kˣ) (hQ : 0 < vAdd v Q) (q : R) (hq : q ∈ I)
    (hQq : algebraMap R K q = (Q : K))
    {a w a' w' : R} {u u' : Kˣ}
    (hu : algebraMap R K a = (u : K)) (hu' : algebraMap R K a' = (u' : K))
    (haw : a * w = q) (haw' : a' * w' = q)
    (h0 : 0 ≤ vAdd v u) (h1 : vAdd v u < vAdd v Q)
    (h0' : 0 ≤ vAdd v u') (h1' : vAdd v u' < vAdd v Q)
    (hcls : QuotientGroup.mk (s := Subgroup.zpowers Q) u = QuotientGroup.mk u')
    (hw1 : IsUnit (1 - w)) (hw1' : IsUnit (1 - w'))
    (hne : algebraMap R K (1 - a) ≠ 0) (hne' : algebraMap R K (1 - a') ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0) :
    tatePtPair (K := K) a w q hq haw hw1 hne hΔ
      = tatePtPair (K := K) a' w' q hq haw' hw1' hne' hΔ := by
  obtain ⟨rfl, rfl⟩ :=
    pair_eq_of_same_class' hinj v Q hQ q hQq hu hu' haw haw' h0 h1 h0' h1' hcls
  rfl

end Class

/-! ## ★出典の紐付け(`.src`) -/

def tatePtPair_eq_of_same_class'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——点は類だけで決まる)",
    sectionId := "genell-def-3-3" }

def normRep_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——単位類の代表元は 1)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
