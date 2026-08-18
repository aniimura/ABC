import ABC3.Found.Arakelov.PicCoverIso

/-!
# Arakelov (B1) 第 103 ブロック —— **生成元による乗法は全単射**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★局所自明性の同型を「乗法」として書く

第 93 ブロックの `unitHomOfSection s` の `app` は **`c ↦ c · s`** である。
★★これが同型であることを言うには、「`s` が生成元なら乗法は全単射」が要る。

★★★**`P ≃ₗ[A] A` があれば、`s := e.symm 1` について `c ↦ c • s` は `e.symm` そのもの**
——`e.symm c = e.symm (c • 1) = c • e.symm 1` だからである。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `generatorOf` | ★生成元 `e.symm 1` |
| `smul_generator_eq` | ★乗法は `e.symm` そのもの |
| `bijective_smul_generator` | ★★★★**生成元による乗法は全単射** |

## ★★★次

★第 100 ブロック(基底の全体で `M_{t·g} ≅ R_{t·g}`)と合わせると、
基底の上で `c ↦ c · s` が全単射になる。
★★あとは第 102 ブロック(被覆で全単射なら同型)を当てれば
`(tilde M)|_{D g} ≅ 𝟙_` が出る。
-/

universe u

namespace ABC3.Found.Arakelov

variable {A : Type u} [CommRing A] {P : Type u} [AddCommGroup P] [Module A P]

/-- ★**生成元**——`P ≃ₗ[A] A` の逆で `1` を送ったもの。 -/
def generatorOf (e : P ≃ₗ[A] A) : P := e.symm 1

/-- ★**乗法は `e.symm` そのものである**。 -/
theorem smul_generator_eq (e : P ≃ₗ[A] A) :
    (fun c : A => c • generatorOf e) = e.symm := by
  funext c
  rw [generatorOf, ← LinearEquiv.map_smul, smul_eq_mul, mul_one]

/-- ★★★★**生成元による乗法は全単射である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが局所自明性の同型の中身である。 -/
theorem bijective_smul_generator (e : P ≃ₗ[A] A) :
    Function.Bijective (fun c : A => c • generatorOf e) := by
  rw [smul_generator_eq]
  exact e.symm.bijective

/-! ## ★出典の紐付け(`.src`) -/

def bijective_smul_generator.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——生成元による乗法は全単射)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
