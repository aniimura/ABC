import ABC3.Found.Arakelov.PicResTrans

/-!
# Arakelov (B1) 第 58 ブロック —— **被覆は引き戻しで輸送される**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★局所自明性の輸送の最後の 1 点

`Opens` の Grothendieck 位相は**点ごと**である:

    S ∈ J U  ⟺  ∀ x ∈ U, ∃ V, ∃ (V ⟶ U), S それ ∧ x ∈ V

★★したがって `Y` の被覆篩 `S` から `X` の開集合 `U` の上の篩

    { V ⟶ U | ∃ W, S (W ⟶ ⊤) ∧ V ≤ f⁻¹W }

を作れば、それは**被覆篩**になる——`x ∈ U` に対し `f(x) ∈ ⊤` なので
`S` から `W ∋ f(x)` が取れ、`x ∈ f⁻¹W ⊓ U` ✅。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `transportSieve` | ★`Y` の篩を `X` の開集合の上の篩へ運ぶ |
| `transportSieve_mem` | ★★★**輸送した篩は被覆篩である** |

## ★★★これで局所自明性の輸送が完成する

    第 58(被覆の輸送)+ 第 54(Beck–Chevalley)+ 第 56(`(f|)^* 𝒪 ≅ 𝒪`)
      + 第 57(制限の推移律)
      ⟹ **IsLocallyTrivial Y L ⟹ IsLocallyTrivial X (f^*_pre L)**
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-- ★**`Y` の被覆篩から `X` の開集合 `U` の上の篩を作る**。 -/
def transportSieve (S : Sieve (⊤ : Y.Opens)) (U : X.Opens) : Sieve U where
  arrows V _ := ∃ W : Y.Opens, S.arrows (homOfLE le_top : W ⟶ ⊤)
    ∧ V ≤ (Opens.map f.base).obj W
  downward_closed := by
    rintro V V' i ⟨W, hW, hV⟩ j
    exact ⟨W, hW, le_trans (leOfHom j) hV⟩

/-- ★★★**輸送した篩は被覆篩である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★`Opens` の Grothendieck 位相は**点ごと**なので、
`x ∈ U` に対し `f(x)` を覆う `W` を取れば `x ∈ f⁻¹W ⊓ U` である。

★★★これが局所自明性の輸送の最後の 1 点である。 -/
theorem transportSieve_mem (S : Sieve (⊤ : Y.Opens))
    (hS : S ∈ (Opens.grothendieckTopology Y) ⊤) (U : X.Opens) :
    transportSieve f S U ∈ (Opens.grothendieckTopology X) U := by
  intro x hx
  obtain ⟨W, g, hg, hW⟩ := hS (f.base x) trivial
  refine ⟨(Opens.map f.base).obj W ⊓ U, homOfLE inf_le_right, ⟨W, ?_, inf_le_left⟩, ?_, hx⟩
  · exact (Subsingleton.elim g (homOfLE le_top : W ⟶ ⊤)) ▸ hg
  · exact hW

/-! ## ★出典の紐付け(`.src`) -/

def transportSieve_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——被覆が引き戻しで輸送されること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
