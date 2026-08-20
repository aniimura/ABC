import ABC3.Found.FrdI.Example43
import Mathlib.CategoryTheory.SingleObj
import Mathlib.CategoryTheory.Localization.Predicate

/-!
# `𝒟` は `𝒞` の base-isomorphism による局所化では**ない**

★★★**本ファイルは反証である**。

`Theorem 3.4, (v)` の `Ψ_Base : 𝒟₁ → 𝒟₂` を組むとき、私は
「`𝒟` は `𝒞` を base-isomorphism で局所化したものである」という**枠**を採り、
`psiBaseLoc` を `[P₁.proj.IsLocalization (baseIsoProp P₁)]` の下で定義した。

★**この仮定は偽である。**

## ★原典は局所化とは言っていない(物理 p.68 を実測)

> which may be composed with a quasi-inverse to the natural equivalence Di

原典は `𝒟_i` を、圏 `(P_i)_{A_𝒟}` の 2-圏の**粗化** `R_i` として復元する
(`Q_i ⊂ R_i`、`R_i ≌ 𝒟_i`)。★**局所化 `𝒞[W^{-1}]` は一度も現れない。**

## ★★反例 —— `Example 4.3` の `𝒞` がそのまま効く

`Example 4.3` では

* `𝒟 = Discrete PUnit`(1 対象・1 射)
* したがって **`𝒞` の射はすべて base-isomorphism**(`ex43_isBaseIso`)

なので `W = baseIsoProp` は**全射類**である。もし `proj : 𝒞 ⥤ 𝒟` が
`W` による局所化なら、`𝒞` の**亜群化が自明圏**でなければならない。

★★しかし `𝒞` の対象 `0 ∈ ℚ` の自己射は `{d ∈ ℕ≥1}` 全体であり、
**Frobenius 次数は亜群化しても死なない** —— 2 進付値
`ℕ≥1 →* Multiplicative ℤ` が `𝒞` から**亜群への関手**を与え、
それが次数 2 の自己射を `1` でない元へ送るからである。

★★★**したがって `𝒪^▷`・`𝒪^×` や Frobenius 次数といった「底が見ない情報」が
局所化には残ってしまい、`𝒟` には落ちない。** これが反証の中身である。

## ★この反証が壊すもの・壊さないもの

| 対象 | 影響 |
|---|---|
| `psiBaseLoc` / `psiBaseLocSquare` / `psiBaseLocUniq` | ★**空虚**(仮定を満たす `𝒞` が無い) |
| `thm_3_4_v_rigid` | ★同上 |
| `proj_inverts_baseIso` | ★**無傷**(局所化の第 1 条件は真) |
| `baseIso_leftFraction_ext` | ★**無傷**(`𝒞` が totally epimorphic なので真) |
| `Theorem 3.4, (v)` 本体 | ★**未証明に戻る** —— 原典の粗化の筋で組み直す |
-/

namespace ABC3.Check.FrdI

open CategoryTheory ABC3.Found.FrdI

/-! ## ★1. 2 進付値 —— `ℕ≥1` から群への単系準同型 -/

/-- ★★`ℕ≥1 →* Multiplicative ℤ` —— `d ↦ ord_2(d)`。
★**行き先が群である**ことが要点で、これで亜群への関手が立つ。 -/
def val2 : ℕ+ →* Multiplicative ℤ where
  toFun d := Multiplicative.ofAdd (((d : ℕ).factorization 2 : ℕ) : ℤ)
  map_one' := by
    show Multiplicative.ofAdd (((1 : ℕ).factorization 2 : ℕ) : ℤ) = 1
    rw [Nat.factorization_one]
    rfl
  map_mul' a b := by
    show Multiplicative.ofAdd ((((a * b : ℕ+) : ℕ).factorization 2 : ℕ) : ℤ)
      = Multiplicative.ofAdd (((a : ℕ).factorization 2 : ℕ) : ℤ)
        * Multiplicative.ofAdd (((b : ℕ).factorization 2 : ℕ) : ℤ)
    have h : ((a * b : ℕ+) : ℕ) = (a : ℕ) * (b : ℕ) := rfl
    rw [h, Nat.factorization_mul a.ne_zero b.ne_zero]
    simp only [Finsupp.add_apply, Nat.cast_add]
    exact ofAdd_add _ _

theorem val2_two : val2 2 = Multiplicative.ofAdd (1 : ℤ) := by
  show Multiplicative.ofAdd ((((2 : ℕ+) : ℕ).factorization 2 : ℕ) : ℤ) = _
  norm_num [Nat.Prime.factorization_self Nat.prime_two]

theorem val2_two_ne_one : val2 2 ≠ 1 := by
  rw [val2_two]
  intro h
  have : (1 : ℤ) = 0 := congrArg Multiplicative.toAdd h
  exact one_ne_zero this

/-! ## ★2. `𝒞 = Ex43` から 1 対象亜群への関手

★`SingleObj M` の合成は `comp x y = y * x`(mathlib)であり、
`Ex43` の合成も `(f ≫ g).deg = g.deg * f.deg` なので**向きが合う**。 -/

/-- ★★★**`Example 4.3` の `𝒞` から 1 対象亜群への関手** —— 射を 2 進付値へ送る。 -/
def ex43Val2 : Ex43 ⥤ SingleObj (Multiplicative ℤ) where
  obj _ := SingleObj.star _
  map f := val2 f.deg
  map_id a := by
    show val2 (𝟙 a : Ex43Hom a a).deg = (1 : Multiplicative ℤ)
    rw [Ex43.id_deg]
    exact val2.map_one
  map_comp {a b c} f g := by
    show val2 (f ≫ g : Ex43Hom a c).deg = val2 g.deg * val2 f.deg
    rw [Ex43.comp_deg]
    exact val2.map_mul _ _

/-- ★`ex43Val2` は**すべての射を同型へ送る** —— 行き先が亜群だから。 -/
theorem ex43Val2_inverts :
    (baseIsoProp ex43P).IsInvertedBy ex43Val2 := fun _ _ _ _ => inferInstance

/-! ## ★3. 反証の本体 -/

/-- ★★★**`proj : 𝒞 ⥤ 𝒟` は base-isomorphism による局所化では**ない****。

★`Example 4.3` の `𝒞` が反例になる。 -/
theorem ex43_not_isLocalization :
    ¬ (ex43P.proj.IsLocalization (baseIsoProp ex43P)) := by
  intro hloc
  haveI := hloc
  -- 局所化なら `ex43Val2` は `𝒟` を経由して分解する
  let G : Discrete PUnit ⥤ SingleObj (Multiplicative ℤ) :=
    Localization.lift ex43Val2 ex43Val2_inverts ex43P.proj
  have hfac : ex43P.proj ⋙ G ≅ ex43Val2 :=
    Localization.fac ex43Val2 ex43Val2_inverts ex43P.proj
  -- 対象 `0 ∈ ℚ` と、その次数 2 の自己射
  let a : Ex43 := Ex43.mk 0
  have hcond : ((((2 : ℕ+) : ℕ) : ℚ)) * a.val ≤ a.val := by
    show ((2 : ℕ) : ℚ) * (0 : ℚ) ≤ (0 : ℚ)
    norm_num
  let f2 : a ⟶ a := ⟨2, hcond⟩
  -- `𝒟` の側では `f2` と `𝟙 a` は同じ射になる
  have hproj : (ex43P.proj ⋙ G).map f2 = (ex43P.proj ⋙ G).map (𝟙 a) := by
    have : ex43P.proj.map f2 = ex43P.proj.map (𝟙 a) := Subsingleton.elim _ _
    exact congrArg G.map this
  -- 自然性から `ex43Val2.map f2 = 𝟙`
  have hid : (ex43P.proj ⋙ G).map f2 = 𝟙 _ := by
    rw [hproj]; exact CategoryTheory.Functor.map_id _ _
  have hnat := hfac.hom.naturality f2
  rw [hid, Category.id_comp] at hnat
  haveI : IsIso (hfac.hom.app a) := inferInstance
  have h2 : hfac.hom.app a ≫ ex43Val2.map f2 = hfac.hom.app a ≫ 𝟙 _ := by
    rw [Category.comp_id]; exact hnat.symm
  have hkey : ex43Val2.map f2 = 𝟙 _ := (cancel_epi (hfac.hom.app a)).mp h2
  -- ところが `ex43Val2.map f2 = val2 2 ≠ 1`
  exact val2_two_ne_one hkey

/-! ## ★4. 局所化の**第 1 条件**は無傷である

★`proj` がすべての base-isomorphism を同型へ送ること(`proj_inverts_baseIso`)は
反証の影響を受けない —— 偽なのは**第 2 条件**(`Construction.lift` が圏同値)である。 -/

/-- ★★`Example 4.3` でも第 1 条件は成り立つ。 -/
theorem ex43_proj_inverts :
    (baseIsoProp ex43P).IsInvertedBy ex43P.proj := fun _ _ _ hf => hf

end ABC3.Check.FrdI
