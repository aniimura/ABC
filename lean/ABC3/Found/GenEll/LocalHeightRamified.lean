import ABC3.Meta.Claim
import Mathlib.NumberTheory.RamificationInertia.Valuation

/-!
# GenEll 第 359 ブロック —— **★★★★★★★`Remark 3.3.1` の実質**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.16。

原文 (GenEll p.16):
> that this definition is independent of the choice of L].

## ★★★★★★★★何が `Remark 3.3.1` の中身なのか

原文は 1 文である——「`E_K` が**潜在的**乗法還元しか持たない場合でも、
乗法還元を持つ有限拡大 `L/K` を取り、`E_K ×_K L` の局所高さを
`L/K` の**分岐指数で割れば** `ℚ` の元として定義できる。
★one verifies immediately that this definition is **independent of the choice of L**」。

★★**「immediately」が畳んでいるのは次の事実である**:

> **`ord_w(x) = e(w/v)·ord_v(x)`**(`x ∈ Kˣ`、`w` は `v` の上にある素点)

★★★これがあれば `ord_w(q)/e(w/v) = ord_v(q)` となり、**右辺は `L` を含まない**ので
`L` の取り方に依らないことが直ちに出る。

## ★★★★★2026-08-26 の在庫調査——**mathlib にあった**

| 段 | mathlib |
|---|---|
| **`v.valuation K x ^ e = w.valuation L (algebraMap K L x)`** | ✅ `IsDedekindDomain.HeightOneSpectrum.valuation_liesOver` |
| `e ≠ 0`(`w` が `v` の上にあるとき) | ✅ `ramificationIdx_ne_zero_of_liesOver` |
| 加法的な位数 `ord_v` | ★無い(`valuation` は乗法的な `WithZero (Multiplicative ℤ)` 値) |

★`Found/GaloisRep/NeronValuation.lean` の `valAdd` は**数体の整数環に特化**しており、
かつ第 320 の `vAdd_algebraMap_eq_valAdd` は**不分岐**を仮定している
——`Remark 3.3.1` はまさに**分岐する場合**なので使えない。
★★そこで一般の Dedekind 環に対する `ordAt` を建て、分岐指数倍の法則を出す。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `ordAt` | ★★加法的な位数 `ord_v(x) ∈ ℤ` |
| `valuation_eq_ofAdd_neg_ordAt` | ★乗法的な付値との対応 |
| `ofAdd_pow` | ★`(ofAdd a)^n = ofAdd (n·a)` |
| `ramificationIdx_ne_zero` | ★`e ≠ 0` |
| `ordAt_liesOver` | ★★★★★★**`ord_w(x) = e·ord_v(x)`** |
| `ordAt_div_ramificationIdx` | ★★★★★★★**`ord_w(x)/e = ord_v(x)`** |
| `ordAt_div_ramificationIdx_indep` | ★★★★★★★★**`Remark 3.3.1`** |
-/

namespace ABC3.Found.GenEll

open IsDedekindDomain IsDedekindDomain.HeightOneSpectrum Ideal.IsDedekindDomain

/-! ## ★★加法的な位数 -/

/-- ★★**素点 `v` における加法的な位数** `ord_v(x) ∈ ℤ`。

★mathlib の `valuation` は `WithZero (Multiplicative ℤ)` に値を取るので、
`x` が単元であることを使って指数を取り出す。 -/
noncomputable def ordAt {A K : Type*} [CommRing A] [IsDedekindDomain A] [Field K] [Algebra A K]
    [IsFractionRing A K] (v : HeightOneSpectrum A) (x : Kˣ) : ℤ :=
  -Multiplicative.toAdd (WithZero.unzero (Valuation.ne_zero_of_unit (v.valuation K) x))

theorem valuation_eq_ofAdd_neg_ordAt {A K : Type*} [CommRing A] [IsDedekindDomain A] [Field K]
    [Algebra A K] [IsFractionRing A K] (v : HeightOneSpectrum A) (x : Kˣ) :
    (v.valuation K) (x : K) = (Multiplicative.ofAdd (-(ordAt v x)) : Multiplicative ℤ) := by
  rw [ordAt, neg_neg]
  exact (WithZero.coe_unzero _).symm

theorem ofAdd_pow (a : ℤ) (n : ℕ) :
    (Multiplicative.ofAdd a) ^ n = Multiplicative.ofAdd ((n : ℤ) * a) := by
  induction n with
  | zero => simp
  | succ k ih =>
      rw [pow_succ, ih]
      show Multiplicative.ofAdd _ * Multiplicative.ofAdd _ = _
      rw [← ofAdd_add]
      push_cast
      ring_nf

theorem ramificationIdx_ne_zero {A : Type*} {B : Type*}
    [CommRing A] [IsDedekindDomain A] [CommRing B] [IsDedekindDomain B] [Algebra A B]
    [Module.IsTorsionFree A B]
    (v : HeightOneSpectrum A) (w : HeightOneSpectrum B) [w.asIdeal.LiesOver v.asIdeal] :
    v.asIdeal.ramificationIdx w.asIdeal ≠ 0 :=
  ramificationIdx_ne_zero_of_liesOver w.asIdeal v.ne_bot

/-! ## ★★★★★★分岐指数倍の法則 -/

/-- ★★★★★★**`ord_w(x) = e(w/v)·ord_v(x)`**(`x ∈ Kˣ`)。

原文 (GenEll p.16):
> that this definition is independent of the choice of L].

★mathlib の `valuation_liesOver`(乗法的な形)を加法的に読み替えただけである。 -/
theorem ordAt_liesOver {A K : Type*} (L : Type*) {B : Type*}
    [CommRing A] [IsDedekindDomain A] [CommRing B] [IsDedekindDomain B] [Algebra A B]
    [Module.IsTorsionFree A B] [Field K] [Field L] [Algebra K L]
    [Algebra A K] [IsFractionRing A K] [Algebra A L] [IsScalarTower A K L]
    [Algebra B L] [IsFractionRing B L] [IsScalarTower A B L]
    (v : HeightOneSpectrum A) (w : HeightOneSpectrum B) [w.asIdeal.LiesOver v.asIdeal]
    (x : Kˣ) :
    ordAt w (Units.map (algebraMap K L).toMonoidHom x)
      = (v.asIdeal.ramificationIdx w.asIdeal : ℤ) * ordAt v x := by
  have h := valuation_liesOver L v w (x : K)
  rw [valuation_eq_ofAdd_neg_ordAt v x] at h
  have hy : ((Units.map (algebraMap K L).toMonoidHom x : Lˣ) : L) = algebraMap K L (x : K) := rfl
  rw [← hy, valuation_eq_ofAdd_neg_ordAt w] at h
  have h2 : ((Multiplicative.ofAdd (-(ordAt v x)) : Multiplicative ℤ)
      : WithZero (Multiplicative ℤ)) ^ (v.asIdeal.ramificationIdx w.asIdeal)
      = ((Multiplicative.ofAdd
          (-((v.asIdeal.ramificationIdx w.asIdeal : ℤ) * ordAt v x)) : Multiplicative ℤ)
          : WithZero (Multiplicative ℤ)) := by
    rw [← WithZero.coe_pow, ofAdd_pow]
    congr 2
    ring
  rw [h2] at h
  have h3 := WithZero.coe_inj.1 h
  have h4 := Multiplicative.ofAdd.injective h3
  linarith

/-- ★★★★★★★**`ord_w(x)/e(w/v) = ord_v(x)`**——右辺は `L` を含まない。 -/
theorem ordAt_div_ramificationIdx {A K : Type*} (L : Type*) {B : Type*}
    [CommRing A] [IsDedekindDomain A] [CommRing B] [IsDedekindDomain B] [Algebra A B]
    [Module.IsTorsionFree A B] [Field K] [Field L] [Algebra K L]
    [Algebra A K] [IsFractionRing A K] [Algebra A L] [IsScalarTower A K L]
    [Algebra B L] [IsFractionRing B L] [IsScalarTower A B L]
    (v : HeightOneSpectrum A) (w : HeightOneSpectrum B) [w.asIdeal.LiesOver v.asIdeal]
    (x : Kˣ) :
    ((ordAt w (Units.map (algebraMap K L).toMonoidHom x) : ℚ))
        / ((v.asIdeal.ramificationIdx w.asIdeal : ℕ) : ℚ)
      = (ordAt v x : ℚ) := by
  have he : ((v.asIdeal.ramificationIdx w.asIdeal : ℕ) : ℚ) ≠ 0 := by
    exact_mod_cast ramificationIdx_ne_zero v w
  rw [ordAt_liesOver L v w x]
  push_cast
  field_simp

/-- ★★★★★★★★**`Remark 3.3.1`**——`ord_w(q)/e` は拡大 `L/K` の取り方に依らない。

原文 (GenEll p.16):
> that this definition is independent of the choice of L].

★2 つの拡大の値がどちらも `ord_v(q)` に等しいので、直ちに一致する。 -/
theorem ordAt_div_ramificationIdx_indep {A K : Type*} (L L' : Type*) {B B' : Type*}
    [CommRing A] [IsDedekindDomain A] [Field K] [Algebra A K] [IsFractionRing A K]
    [CommRing B] [IsDedekindDomain B] [Algebra A B] [Module.IsTorsionFree A B]
    [Field L] [Algebra K L] [Algebra A L] [IsScalarTower A K L]
    [Algebra B L] [IsFractionRing B L] [IsScalarTower A B L]
    [CommRing B'] [IsDedekindDomain B'] [Algebra A B'] [Module.IsTorsionFree A B']
    [Field L'] [Algebra K L'] [Algebra A L'] [IsScalarTower A K L']
    [Algebra B' L'] [IsFractionRing B' L'] [IsScalarTower A B' L']
    (v : HeightOneSpectrum A) (w : HeightOneSpectrum B) [w.asIdeal.LiesOver v.asIdeal]
    (w' : HeightOneSpectrum B') [w'.asIdeal.LiesOver v.asIdeal]
    (x : Kˣ) :
    ((ordAt w (Units.map (algebraMap K L).toMonoidHom x) : ℚ))
        / ((v.asIdeal.ramificationIdx w.asIdeal : ℕ) : ℚ)
      = ((ordAt w' (Units.map (algebraMap K L').toMonoidHom x) : ℚ))
        / ((v.asIdeal.ramificationIdx w'.asIdeal : ℕ) : ℚ) := by
  rw [ordAt_div_ramificationIdx L v w x, ordAt_div_ramificationIdx L' v w' x]

/-! ## ★出典の紐付け(`.src`) -/

def ordAt_liesOver.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 16, item := "Remark 3.3.1(局所高さの ℚ への拡張)",
    sectionId := "genell-rem-3-3-1" }

def ordAt_div_ramificationIdx_indep.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 16, item := "Remark 3.3.1(L の取り方に依らないこと)",
    sectionId := "genell-rem-3-3-1" }

end ABC3.Found.GenEll
