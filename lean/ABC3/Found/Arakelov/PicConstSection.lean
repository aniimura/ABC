import ABC3.Found.Arakelov.PicTensorUniq

/-!
# Arakelov (B1) 第 87 ブロック —— **定数切断のテンソルは定数切断**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★比較射を同定する核心の恒等式

    tensorSection (const m) (const n)  =  const (m ⊗ₜ n)

★★これが「比較射は `M ⊗_R N` の局所化の標準の写像である」ことの中身である。

## ★★機構 —— 点ごとに第 80 ブロック

`toOpenₗ R M U m` は**定数切断** `x ↦ mk m 1` である(`rfl` で確認)。
★したがって点ごとに

    localizedTensorEquiv (mk m 1 ⊗ₜ mk n 1) = mk (m ⊗ₜ n) 1

——第 80 ブロックの `localizedTensorEquiv_mk_one` そのものである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `toOpenL_apply` | ★定数切断の値は `mk m 1`(`rfl`) |
| `tensorSection_toOpenL` | ★★★★**定数切断のテンソルは定数切断** |

## ★★★残り

★この恒等式と第 86 ブロックの `tensorLocalization_ext` を合わせると、
比較射の `D(f)` 成分が**局所化の標準写像**と一致し、したがって同型であることが出る。
★★あとは茎(第 77・78)で結論するだけである。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace TensorProduct StructureSheaf

variable (R : Type u) [CommRing R] (M N : Type u)
  [AddCommGroup M] [Module R M] [AddCommGroup N] [Module R N]
  (U : Opens (PrimeSpectrum.Top R))

/-- ★**定数切断の値は `mk m 1` である**(`rfl`)。 -/
theorem toOpenL_apply (m : M) (x : U) : (toOpenₗ R M U m).1 x = LocalizedModule.mk m 1 := rfl

/-- ★★★★**定数切断のテンソルは定数切断である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが比較射を「局所化の標準写像」と同定する核心である
——点ごとに第 80 ブロックの `localizedTensorEquiv_mk_one` を当てるだけ。 -/
theorem tensorSection_toOpenL (m : M) (n : N) :
    tensorSection R M N U (toOpenₗ R M U m) (toOpenₗ R N U n)
      = toOpenₗ R (M ⊗[R] N) U (m ⊗ₜ[R] n) := by
  apply Subtype.ext
  funext x
  exact localizedTensorEquiv_mk_one R _ M N m n

/-! ## ★出典の紐付け(`.src`) -/

def tensorSection_toOpenL.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——定数切断のテンソルは定数切断)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
