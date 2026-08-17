import ABC3.Meta.Claim
import Mathlib.AlgebraicGeometry.IdealSheaf.Basic

/-!
# [GenEll] Definition 1.1 の足場 —— イデアル層の**茎**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> There is an evident notion of tensor product of arithmetic line bundles on X. The

## ★★★なぜ茎で定義し直すのか —— 摩擦の除去

有効 Cartier 因子を **「各点のまわりのアフィン開集合の上で非零因子 1 つで生成」**
と定義すると(`CartierPullback.lean` の `IsEffectiveCartier`)、
「積で閉じる」を示すときに **2 つの摩擦**が出る:

1. **`X.affineOpens` と `X.Opens` の型階層** —— `homOfLE h` は `h` の型で
   行き先の圏が変わる
2. **`X.basicOpen a = X.basicOpen b` に沿った輸送** —— 等号が開集合の間なので
   `Γ(X, ·)` の型が変わり、`rw` が motive で失敗する

★★**茎で定義すると、どちらも消える。**
茎は点だけに依り、開集合を名指ししないからである。
`(I·J)_x = I_x · J_x` は `Ideal.map_mul` 1 本で出る。

★これで「正面から要ると思ったものが要らなかった」の **5 例目**になる
(`Cartier` / scheme 論的像 / NCBelyi の逐語 / 複素解析空間 に続く)。

## ★mathlib 実測(2026-08-17)

`IdealSheafData` に**茎は無い**(`IdealSheaf/` 配下を `stalk` で検索して、
出てくるのは `Subscheme.lean` の `stalkFunctor` の使用のみ)。

## ★選択を使わない定義

各アフィン開集合ごとの像の**上限**として定める。
★これなら well-defined 性が定義に組み込まれ、証明が要らない。
そのうえで「どのアフィン開集合で計算しても同じ」(`stalk_eq_map`)を示す。
-/

namespace AlgebraicGeometry.Scheme.IdealSheafData

open CategoryTheory TopologicalSpace

variable {X : Scheme} (I J : X.IdealSheafData) (x : X)

/-- ★**イデアル層の茎** —— 各アフィン開集合での像の上限。

★選択を使わないので well-defined 性の証明が要らない。 -/
noncomputable def stalk : Ideal (X.presheaf.stalk x) :=
  ⨆ U : {U : X.affineOpens // x ∈ U.1},
    (I.ideal U.1).map (X.presheaf.germ U.1.1 x U.2).hom

variable {I J x}

/-- ★各アフィン開集合の像は茎に含まれる。 -/
theorem le_stalk {U : X.affineOpens} (hx : x ∈ U.1) :
    (I.ideal U).map (X.presheaf.germ U.1 x hx).hom ≤ I.stalk x :=
  le_iSup (fun U : {U : X.affineOpens // x ∈ U.1} =>
    (I.ideal U.1).map (X.presheaf.germ U.1.1 x U.2).hom) ⟨U, hx⟩

/-- ★★**茎はどのアフィン開集合で計算しても同じ**。

★証明の要は「2 つのアフィン開集合の共通部分に、`x` を含むアフィン開集合が取れる」
(`Scheme.isBasis_affineOpens`)ことと、
`IdealSheafData.map_ideal`(制限)+ `germ_res`(茎への写像は制限と両立)である。 -/
theorem stalk_eq_map {U : X.affineOpens} (hx : x ∈ U.1) :
    I.stalk x = (I.ideal U).map (X.presheaf.germ U.1 x hx).hom := by
  refine le_antisymm ?_ (le_stalk hx)
  refine iSup_le fun V => ?_
  -- `U ⊓ V` の中に `x` を含むアフィン開集合 `W` を取る
  obtain ⟨_, ⟨W, hW, rfl⟩, hxW, hWUV⟩ :=
    X.isBasis_affineOpens.exists_subset_of_mem_open
      (show x ∈ (U.1 ⊓ V.1.1 : X.Opens) from ⟨hx, V.2⟩) (U.1 ⊓ V.1.1).2
  have hWU : (⟨W, hW⟩ : X.affineOpens) ≤ U := fun z hz => (hWUV hz).1
  have hWV : (⟨W, hW⟩ : X.affineOpens) ≤ V.1 := fun z hz => (hWUV hz).2
  -- どちらの側も `W` での像に一致する
  have key : ∀ (T : X.affineOpens) (h : (⟨W, hW⟩ : X.affineOpens) ≤ T) (hxT : x ∈ T.1),
      (I.ideal T).map (X.presheaf.germ T.1 x hxT).hom
        = (I.ideal ⟨W, hW⟩).map (X.presheaf.germ (⟨W, hW⟩ : X.affineOpens).1 x hxW).hom := by
    intro T h hxT
    rw [← I.map_ideal h, Ideal.map_map, ← CommRingCat.hom_comp]
    erw [X.presheaf.germ_res]
  rw [key V.1 hWV V.2, ← key U hWU hx]

/-- ★★**茎は積を保つ** —— `(I·J)_x = I_x · J_x`。

★★開被覆版ではここに 2 つの摩擦があった。**茎では `Ideal.map_mul` 1 本である。** -/
theorem stalk_mul : (I * J).stalk x = I.stalk x * J.stalk x := by
  obtain ⟨_, ⟨U, hU, rfl⟩, hxU, -⟩ :=
    X.isBasis_affineOpens.exists_subset_of_mem_open (Set.mem_univ x) isOpen_univ
  rw [stalk_eq_map (U := ⟨U, hU⟩) hxU, stalk_eq_map (U := ⟨U, hU⟩) hxU,
    stalk_eq_map (U := ⟨U, hU⟩) hxU]
  show ((I.ideal ⟨U, hU⟩) * (J.ideal ⟨U, hU⟩)).map _ = _
  rw [Ideal.map_mul]

/-- ★**空因子の茎は全体**。 -/
@[simp] theorem stalk_top : (⊤ : X.IdealSheafData).stalk x = ⊤ := by
  obtain ⟨_, ⟨U, hU, rfl⟩, hxU, -⟩ :=
    X.isBasis_affineOpens.exists_subset_of_mem_open (Set.mem_univ x) isOpen_univ
  rw [stalk_eq_map (U := ⟨U, hU⟩) hxU]
  show (⊤ : Ideal Γ(X, U)).map _ = ⊤
  rw [Ideal.map_top]

end AlgebraicGeometry.Scheme.IdealSheafData

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory

/-! ## ★★茎による有効 Cartier 因子 -/

/-- ★★**有効 Cartier 因子(茎版)** —— 各点の茎が**非零因子 1 つ**で生成される。

★Stacks 01WQ はこの形と開被覆の形の両方を挙げている。
★★**積で閉じることが 3 行で出る**のがこの形の値打ちである。 -/
def IsEffectiveCartierStalk {X : Scheme} (I : X.IdealSheafData) : Prop :=
  ∀ x : X, ∃ f : X.presheaf.stalk x,
    I.stalk x = Ideal.span {f} ∧ f ∈ nonZeroDivisors (X.presheaf.stalk x)

/-- ★★★**有効 Cartier 因子は積で閉じる**(= 直線束のテンソル積の因子側)。

原文 (GenEll p.3):
> There is an evident notion of tensor product of arithmetic line bundles on X. The

★★**開被覆版では 2 つの摩擦で組み上がらなかったものが、茎では 3 行で出る。** -/
theorem isEffectiveCartierStalk_mul {X : Scheme} {I J : X.IdealSheafData}
    (hI : IsEffectiveCartierStalk I) (hJ : IsEffectiveCartierStalk J) :
    IsEffectiveCartierStalk (I * J) := by
  intro x
  obtain ⟨f, hf, hfn⟩ := hI x
  obtain ⟨g, hg, hgn⟩ := hJ x
  exact ⟨f * g, by
    rw [Scheme.IdealSheafData.stalk_mul, hf, hg, Ideal.span_singleton_mul_span_singleton],
    Submonoid.mul_mem _ hfn hgn⟩

/-- ★**空因子は有効 Cartier**(非空虚性の witness、茎版)。 -/
theorem isEffectiveCartierStalk_top (X : Scheme) :
    IsEffectiveCartierStalk (⊤ : X.IdealSheafData) := by
  intro x
  refine ⟨1, ?_, Submonoid.one_mem _⟩
  rw [Scheme.IdealSheafData.stalk_top, Ideal.span_singleton_one]

/-! ## ★出典の紐付け(`.src`) -/

def IsEffectiveCartierStalk.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1(テンソル積の因子側——モノイド構造のみ)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.GenEll
