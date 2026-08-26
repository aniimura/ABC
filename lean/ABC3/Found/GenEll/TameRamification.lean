import ABC3.Meta.Claim
import Mathlib.RingTheory.LocalRing.ResidueField.Defs
import Mathlib.RingTheory.LocalRing.Basic
import Mathlib.Algebra.CharP.Basic

/-!
# [GenEll] Proposition 1.7 の "elementary claim" —— **馴分岐**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.10。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

## ★★★★★★mathlib に馴分岐は無い(2026-08-26 実測)

REPL で名前を引いた結果:

| 探したもの | mathlib |
|---|---|
| `Algebra.IsUnramifiedAt` | ✅ ある |
| `IsTamelyRamified` / `Algebra.IsTamelyRamified` / `Ideal.IsTamelyRamified` / `IsTame` | ❌ **どれも無い** |

★したがって**建てる**しかない。本ファイルはその最初の一段である。

## ★★★何を取るか —— 第 379 の仮説は馴分岐そのものである

`DifferentKummer.lean` の `mem_differentIdeal_of_isUnit_natCast`(第 379)は
**`IsUnit (n : B)`** を仮定していた。★これが**馴分岐の言い換え**であることを示す:

> 局所環 `B` で `IsUnit (n : B)` ⟺ **剰余標数が `n` を割らない**

★★これで「野生か馴か」が **`n` が単元か**という 1 つの述語に集約される。
★★★原文の段 1(野生と馴に分ける)は、形式化の側では**この 1 行の場合分け**になる。
-/

namespace ABC3.Found.GenEll

open IsLocalRing

/-! ## ★★★★馴分岐の定義 -/

/-- ★★★★**馴分岐の次数**(局所環の言葉で)—— 剰余標数が `n` を割らない。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

★原文が『tamely ramified』と呼ぶものを、**次数の側だけ**取り出した形である
(剰余体の分離性は `Algebra.IsSeparable` として別に持てる)。 -/
def IsTameDegree (B : Type*) [CommRing B] [IsLocalRing B] (n : ℕ) : Prop :=
  ¬ (ringChar (ResidueField B) ∣ n)

/-! ## ★★★★★`IsUnit (n : B)` は馴分岐そのもの -/

/-- ★★★★★**局所環では `IsUnit (n : B)` ⟺ 剰余標数が `n` を割らない**。

★機構は 3 行——単元であることは極大イデアルに入らないこと、
極大イデアルに入ることは剰余体で `0` になること、
剰余体で `0` になることは標数が割ることである。 -/
theorem isUnit_natCast_iff (B : Type*) [CommRing B] [IsLocalRing B] (ell n : ℕ)
    (hchar : ringChar (ResidueField B) = ell) :
    IsUnit (n : B) ↔ ¬ (ell ∣ n) := by
  have hres : ∀ x : B, residue B x = 0 ↔ x ∈ maximalIdeal B := fun x =>
    ⟨fun h => (Ideal.Quotient.eq_zero_iff_mem (I := maximalIdeal B) (a := x)).mp h,
     fun h => (Ideal.Quotient.eq_zero_iff_mem (I := maximalIdeal B) (a := x)).mpr h⟩
  have hcast : ((n : ℕ) : ResidueField B) = residue B (n : B) := by simp
  rw [← notMem_maximalIdeal, ← hres, ← hcast, ringChar.spec, hchar]

/-- ★★★★★★**第 379 の仮説は馴分岐である**(言い換え)。

★これで `DifferentKummer.lean` の `mem_differentIdeal_of_isUnit_natCast` が
「**馴分岐なら `p·O_L ⊆ 𝔡`**」を言っていることが型で見える。 -/
theorem isUnit_natCast_iff_isTameDegree (B : Type*) [CommRing B] [IsLocalRing B] (n : ℕ) :
    IsUnit (n : B) ↔ IsTameDegree B n :=
  isUnit_natCast_iff B _ n rfl

/-- ★不分岐(次数 `1`)は馴分岐である。 -/
theorem isTameDegree_one (B : Type*) [CommRing B] [IsLocalRing B]
    (hchar : ringChar (ResidueField B) ≠ 1) : IsTameDegree B 1 := by
  intro h
  exact hchar (Nat.dvd_one.mp h)

/-! ## ★出典の紐付け(`.src`) -/

def IsTameDegree.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(段 1——馴分岐の定義。次数の側)",
    sectionId := "genell-prop-1-7" }

def isUnit_natCast_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(段 1——IsUnit (n : B) が馴分岐であること)",
    sectionId := "genell-prop-1-7" }

end ABC3.Found.GenEll
