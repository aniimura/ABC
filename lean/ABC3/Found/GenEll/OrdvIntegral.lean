import ABC3.Found.GenEll.ArithDiv

/-!
# [GenEll] §1 —— `ord_v` による整性の判定(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where cv ∈ Z if v ∈ V(F )non and cv ∈ R if v ∈ V(F )arc [cf. [Szp], §1.1]. Here, if

## ★何を取るか

> **`y ∈ 𝓞_F` ⟺ すべての有限素点で `ord_v(y) ≥ 0`**

★`ArithDiv.lean` は `ordv_algebraMap_nonneg`(整元なら `ord_v ≥ 0`)を持っていたが、
**逆向き**を持っていなかった。逆向きは「分母を潰す」議論に要る。

★mathlib の `IsDedekindDomain.HeightOneSpectrum.mem_integers_of_valuation_le_one`
(すべての付値が 1 以下なら整)に、`ord_v ≥ 0 ⟺ v(y) ≤ 1` を繋ぐだけである。

## ★★符号の 3 度目の確認になる

`ordv` の符号は 2026-08-17 に一度訂正している。
★本ファイルの `valuation_le_one_iff_ordv_nonneg` は
「付値が小さい ⟺ `ord` が大きい」という**向き**を固定するので、
符号が逆なら整性の判定が反転して破れる。
-/

namespace ABC3.Found.GenEll

open NumberField IsDedekindDomain

variable {F : Type*} [Field F] [NumberField F]

/-- ★**`v(f) ≤ 1 ⟺ ord_v(f) ≥ 0`**。

★`ordv` は `toAdd` の**符号を反転したもの**なので、
付値が小さいことと `ord` が大きいことが対応する。 -/
theorem valuation_le_one_iff_ordv_nonneg (v : FinitePlace F) (f : Fˣ) :
    v.valuation F ((f : F)) ≤ 1 ↔ 0 ≤ ordv v f := by
  have hne : v.valuation F ((f : F)) ≠ 0 :=
    (Valuation.ne_zero_iff (v.valuation F)).2 (Units.ne_zero f)
  have hcoe : v.valuation F ((f : F))
      = ((WithZero.unzero hne : Multiplicative ℤ) : WithZero (Multiplicative ℤ)) :=
    (WithZero.coe_unzero hne).symm
  rw [hcoe, show (1 : WithZero (Multiplicative ℤ))
      = ((1 : Multiplicative ℤ) : WithZero (Multiplicative ℤ)) from rfl,
    WithZero.coe_le_coe]
  simp only [ordv, neg_nonneg]
  exact ⟨fun h => h, fun h => h⟩

/-- ★★**`ord_v ≥ 0` がすべての有限素点で成り立てば整元である**。

★`ArithDiv.lean` の `ordv_algebraMap_nonneg` の**逆**である。
mathlib の `mem_integers_of_valuation_le_one` に上の同値を繋いだ。 -/
theorem mem_range_algebraMap_of_ordv_nonneg (f : Fˣ)
    (h : ∀ v : FinitePlace F, 0 ≤ ordv v f) :
    (f : F) ∈ (algebraMap (𝓞 F) F).range := by
  refine IsDedekindDomain.HeightOneSpectrum.mem_integers_of_valuation_le_one
    (R := 𝓞 F) (K := F) ((f : F)) (fun v => ?_)
  exact (valuation_le_one_iff_ordv_nonneg v f).2 (h v)

/-! ## ★出典の紐付け(`.src`) -/

def valuation_le_one_iff_ordv_nonneg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(算術因子 ADiv(F))",
    sectionId := "genell-adiv" }

end ABC3.Found.GenEll
