import ABC3.Found.Arakelov.PicSectionEquiv

/-!
# Arakelov (B1) 第 128 ブロック —— **点ごとの自明化から `IsLocallyTrivial`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★局所自明性の定義の側を片付ける

`IsLocallyTrivial` は**すべての開集合 `U`** について被覆篩を要求する。
★★しかし実際に手元にあるのは「**各点に自明化する近傍がある**」という形である。

★★★**両者は同値**である——篩を「自明化する近傍に含まれる開集合」で生成すればよい。
★下方閉性は第 57 ブロック(制限の推移律、`rfl`)が保証する。

## ★★機構

| 段 | 内容 |
|---|---|
| 被覆であること | ★点ごとに近傍を取り `U` と交わる |
| 篩の元で自明 | ★★`restrictOnFunctor` で制限(第 57 の `restrict_trans` / `restrictOnUnit`) |

★**`rw` で 2 つの `rfl` 等式を当てるだけ**で閉じた。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isLocallyTrivial_of_pointwise` | ★★★★★★**点ごとの自明化 ⟹ `IsLocallyTrivial`** |

## ★★★これで残るのは「各点に自明化する近傍がある」ことだけ

`tilde M`(`M` 可逆)については第 76(`D(g)` で `M_g ≅ R_g`)+
第 115(切断から同型)+ 第 127(全単射)で出る。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} (P : X.PresheafOfModules)

/-- ★★★★★★**点ごとの自明化から `IsLocallyTrivial` が出る**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★篩を「自明化する近傍に含まれる開集合」で生成すればよい
——下方閉性は第 57 ブロック(制限の推移律)が保証する。 -/
theorem isLocallyTrivial_of_pointwise
    (h : ∀ x : X, ∃ V : X.Opens, x ∈ V ∧
      Nonempty ((restrictPresheafFunctor X V).obj P ≅ 𝟙_ (PresheafModulesOn X V))) :
    IsLocallyTrivial X P := by
  intro U
  refine ⟨Sieve.generate (fun (W : X.Opens) (_ : W ⟶ U) =>
      ∃ V : X.Opens, W ≤ V ∧ Nonempty ((restrictPresheafFunctor X V).obj P
        ≅ 𝟙_ (PresheafModulesOn X V))), ?_, ?_⟩
  · intro x hx
    obtain ⟨V, hxV, hV⟩ := h x
    refine ⟨V ⊓ U, homOfLE inf_le_right, ?_, ⟨hxV, hx⟩⟩
    exact ⟨V ⊓ U, 𝟙 _, homOfLE inf_le_right, ⟨V, inf_le_left, hV⟩, rfl⟩
  · rintro W i ⟨Z, u, v, ⟨V, hZV, ⟨eV⟩⟩, rfl⟩
    have hWV : W ≤ V := le_trans (leOfHom u) hZV
    refine ⟨?_⟩
    have hiso := (restrictOnFunctor hWV).mapIso eV
    rw [restrict_trans hWV P, restrictOnUnit hWV] at hiso
    exact hiso

/-! ## ★出典の紐付け(`.src`) -/

def isLocallyTrivial_of_pointwise.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——点ごとの自明化から IsLocallyTrivial)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
