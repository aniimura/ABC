/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.AwayRange
import ABC3.Found.GenEll.ChartSurjective
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★段 E3c（完成形）—— 座標の条件だけが残る（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★これは何か

段 E3c は「環準同型 `A⁰_{x_i} → Γ(X, X_{s_i} ⊓ V)` が**全射**」である。
本ファイルはそれを

> **ある有限集合 `T`** があって、`T` の各元 `g` に対し
> `g · s_i = s_j` となる座標 `j` が取れるなら、**全射である**

という形に落とす。★★すなわち残るのは「**座標の選び方**」だけであり、
しかもそれは**有限個の確認**である。

## ★★★★★組み立て —— 4 本の合流

| 段 | 出典 | 役割 |
|---|---|---|
| 有限生成 | `§9-832` `exists_finset_generating` | 試すべき元は有限個（段 E3b） |
| 全射判定 | `§9-833` `surjective_of_adjoin_le_range` | 生成元が像に入れば全射 |
| 比の一意性 | `§9-837` `mem_range_awayToSectionHom_of_mul_eq` | `g·s_i = s_j` なら像に入る |
| 定数 | 本ファイル `awayToSectionHom_awayConst` | ★**係数環の像は自動で入る** |

★4 本目が本ファイルの新しい部分である——`A⁰_{x_i}` の**次数 0 の元**（`C r / x_i^0`）が
`φ r` の制限に写ることを見れば、底環の条件（`§9-833` の `hφ`）が自動で埋まる。

## ★測定の記録

★`awayToSectionHom … (awayConst r) = res (φ r)` の証明は
`homogRatio` を開いて `pow_zero`・`eval₂_C` を使うだけである（3 行）。
★★`appLE` の合成則（`res ∘ f.appLE U V = f.appLE U W`）は mathlib に見当たらなかったので
`appLE_res` を自前で置いた（`Opens` の射が `Subsingleton` なのが効く）。

## ★残っている段（明示）

★★★**座標の選び方**だけである: `s : Fin (N+1) → M(⊤)` を
「各チャートの試験元 `g` に対し `a_g` を含む」ように取る段。
★`a_g` は `§9-831`（段 E3a-3）が与える（ただし `n` 乗の形なので
`M` を `M^{⊗n}`・`s_i` を `s_i^{⊗n}` に取り替える、`§9-825`）。
★★有限個（チャート有限 × 試験元有限）なので `N` を大きく取れば済む。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MonoidalCategory ABC3.Found.Arakelov
open MvPolynomial

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}} {A : CommRingCat.{0}}

/-! ## ★`appLE` は制限と合成できる -/

/-- ★**`appLE` を更に制限したものは、より小さい開への `appLE` である**。

★mathlib に見当たらなかったので自前で置いた。`Opens` の射が `Subsingleton` なのが効く。 -/
theorem appLE_res {Y : Scheme.{0}} (f : X ⟶ Y) (U : Y.Opens) (V W : X.Opens)
    (h : V ≤ f ⁻¹ᵁ U) (hWV : W ≤ V) (x : (Γ(Y, U) : Type)) :
    X.presheaf.map (homOfLE hWV).op (f.appLE U V h x)
      = f.appLE U W (hWV.trans h) x := by
  show X.presheaf.map (homOfLE hWV).op
      ((f.app U ≫ X.presheaf.map (homOfLE h).op) x)
    = (f.app U ≫ X.presheaf.map (homOfLE (hWV.trans h)).op) x
  rw [CommRingCat.comp_apply, CommRingCat.comp_apply,
    show (homOfLE (hWV.trans h) : W ⟶ f ⁻¹ᵁ U) = homOfLE hWV ≫ homOfLE h from
      Subsingleton.elim _ _, op_comp, X.presheaf.map_comp, CommRingCat.comp_apply]

/-! ## ★★★★★係数環の像は自動で入る -/

/-- ★**定数を `A⁰_{x_i}` の元と見る**（次数 `0`、`C r / x_i^0`）。 -/
noncomputable def awayConst (N : ℕ) (R : Type) [CommRing R] (i : Fin (N + 1)) (r : R) :
    HomogeneousLocalization.Away
      (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R) (MvPolynomial.X i) :=
  HomogeneousLocalization.Away.mk _ (MvPolynomial.isHomogeneous_X R i) 0 (MvPolynomial.C r)
    (by simp)

/-- ★★★★**定数の像は `φ r` の制限である**——底環の条件が自動で埋まる。

★機構は `homogRatio` を開いて `pow_zero`・`eval₂_C` を使うだけである。 -/
theorem awayToSectionHom_awayConst (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, V) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1)) (r : R) :
    awayToSectionHom M V e φ s i (awayConst N R i r)
      = X.presheaf.map (homOfLE (inf_le_right :
          nonVanishing M (s i) ⊓ V ≤ V)).op (φ r) := by
  show homogRatio M V e φ s i 0 (MvPolynomial.C r) = _
  rw [homogRatio, pow_zero, mul_one]
  congr 1
  show MvPolynomial.eval₂ φ (fun j => trivValue M V e (s j)) (MvPolynomial.C r) = φ r
  rw [MvPolynomial.eval₂_C]

/-- ★★★★★**底環の像は `A⁰_{x_i}` の像に含まれる**（`φ` が底環の像を覆えば）。 -/
theorem range_appLE_subset (f : X ⟶ Spec A)
    (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, V) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1))
    (hφ : ∀ r : (Γ(Spec A, (⊤ : (Spec A).Opens)) : Type),
      f.appLE ⊤ V (by simp) r ∈ Set.range φ) :
    Set.range (f.appLE ⊤ (nonVanishing M (s i) ⊓ V) (by simp)).hom
      ⊆ Set.range (awayToSectionHom M V e φ s i) := by
  rintro _ ⟨r, rfl⟩
  obtain ⟨r', hr'⟩ := hφ r
  refine ⟨awayConst N R i r', ?_⟩
  rw [awayToSectionHom_awayConst, hr']
  exact appLE_res f ⊤ V _ (by simp) inf_le_right r

/-! ## ★★★★★★★★★★段 E3c（完成形） -/

/-- ★★★★★★★★★★**段 E3c（完成形）** —— 残るのは座標の選び方だけ。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

`f : X ⟶ Spec A` が有限型、`X_{s_i} ⊓ V` がアフィン、`φ` が底環の像を覆うとき、
★**ある有限集合 `T`** があって

> `T` の各元 `g` に対し `g · s_i = s_j` となる座標 `j` が取れる

なら `A⁰_{x_i} → Γ(X, X_{s_i} ⊓ V)` は**全射**である。

★★これで段 E3c が「**有限個の座標条件**」に落ちた。
★★★その条件は `§9-831`（段 E3a-3——大域化）が供給する
——ただし `n` 乗の形なので `M` を `M^{⊗n}`・`s_i` を `s_i^{⊗n}` に取り替える（`§9-825`）。 -/
theorem exists_finset_surjective_of_coords (f : X ⟶ Spec A) [LocallyOfFiniteType f]
    (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, V) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1))
    (haff : IsAffineOpen (nonVanishing M (s i) ⊓ V))
    (hφ : ∀ r : (Γ(Spec A, (⊤ : (Spec A).Opens)) : Type),
      f.appLE ⊤ V (by simp) r ∈ Set.range φ) :
    ∃ T : Finset (Γ(X, nonVanishing M (s i) ⊓ V) : Type),
      (∀ g ∈ T, ∃ j : Fin (N + 1),
        g * X.presheaf.map (homOfLE (inf_le_right :
              nonVanishing M (s i) ⊓ V ≤ V)).op (trivValue M V e (s i))
          = X.presheaf.map (homOfLE (inf_le_right :
              nonVanishing M (s i) ⊓ V ≤ V)).op (trivValue M V e (s j))) →
        Function.Surjective (awayToSectionHom M V e φ s i) := by
  obtain ⟨T, hT⟩ := exists_finset_surjective_criterion f haff
  refine ⟨T, fun hgen => hT _ (range_appLE_subset f M V e φ s i hφ) ?_⟩
  intro g hg
  obtain ⟨j, hj⟩ := hgen g hg
  exact mem_range_awayToSectionHom_of_mul_eq M V e φ s i j g hj

/-! ## ★出典の紐付け(`.src`) -/

def appLE_res.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(appLE は制限と合成できる)",
    sectionId := "genell-prop-1-4" }

def awayToSectionHom_awayConst.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(定数の像は φ r の制限——底環の条件は自動)",
    sectionId := "genell-prop-1-4" }

def exists_finset_surjective_of_coords.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(段 E3c 完成形——残るのは座標の選び方だけ)",
    sectionId := "genell-prop-1-4" }

def exists_finset_surjective_of_coords.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_finset_generating(段 E3b、§9-832)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_finset_generating") 2,
    .citation "[ABC3]" "mem_range_awayToSectionHom_of_mul_eq(段 E3c の橋、§9-837)"
      (.inProject "ABC3" "ABC3.Found.GenEll.mem_range_awayToSectionHom_of_mul_eq") 2,
    .citation "[ABC3]" "exists_global_section_of_localData(段 E3a-3、§9-831)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_global_section_of_localData") 3,
    .implicitStep
      ("★appLE の合成則(res ∘ f.appLE U V = f.appLE U W)は mathlib に見当たらなかったので " ++
       "appLE_res を自前で置いた(Opens の射が Subsingleton なのが効く、2026-08-28 実測)") 3,
    .implicitStep
      ("★★★**残るのは座標の選び方だけ**である: s : Fin (N+1) → M(⊤) を" ++
       "「各チャートの試験元 g に対し a_g を含む」ように取る段。" ++
       "★a_g は §9-831 が与える(n 乗の形なので M を M^{⊗n}・s_i を s_i^{⊗n} に取り替える、§9-825)。" ++
       "★★有限個(チャート有限 × 試験元有限)なので N を大きく取れば済む") 7 ]

end ABC3.Found.GenEll
