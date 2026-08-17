import ABC3.Found.GenEll.ArchPoint

/-!
# Arakelov (C1) の第一段 —— **ℂ-点における切断の評価**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★これは `Interface/Arakelov/ArcSpace.lean` の `evalAffine` である

`ArcSpaceData`(C1)は 3 つを要求する:

| 場 | 内容 | 本ファイル |
|---|---|---|
| `Arc X` | 複素点の集合 | ★`complexPoints X`(既存) |
| `evalAffine` | アフィンでの切断の評価 `a(p) ∈ ℂ` | ★★**本ファイル** |
| `topology` | 各点収束・開埋め込みで誘導 | 次段 |

## ★★なぜ評価が要るのか

★★★**位相を固定するのは評価写像である。**
`ArcSpaceData` の `topology_affine` は

    topology (Spec A) = induced (evalAffine A) Pi.topologicalSpace

と言っている。★これが無いと**離散位相で埋まってしまう**(退化検査で実測)。

## ★★構成

`p : Spec ℂ ⟶ Spec A` に対し、`Γ` を取ると `Γ(Spec A) → Γ(Spec ℂ)`、
すなわち `ΓSpecIso` で `A → ℂ` を得る。★これが評価である。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory

/-! ## ★★アフィンでの評価 -/

/-- ★★**ℂ-点 `p : Spec ℂ ⟶ Spec A` が定める環準同型 `A → ℂ`**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★`Spec` は充満忠実なので、アフィンの間の射は環準同型と 1 対 1 に対応する
(`AlgebraicGeometry.Spec.preimage`)。 -/
noncomputable def evalHom (A : CommRingCat.{0})
    (p : Spec (CommRingCat.of ℂ) ⟶ Spec A) : A ⟶ CommRingCat.of ℂ :=
  Spec.preimage p

/-- ★★**切断 `a : A` の点 `p` における値** `a(p) ∈ ℂ`。 -/
noncomputable def evalAffine (A : CommRingCat.{0})
    (p : Spec (CommRingCat.of ℂ) ⟶ Spec A) (a : A) : ℂ :=
  (evalHom A p).hom a

/-! ## ★評価は環準同型である -/

@[simp] theorem evalAffine_one (A : CommRingCat.{0})
    (p : Spec (CommRingCat.of ℂ) ⟶ Spec A) : evalAffine A p 1 = 1 :=
  map_one (evalHom A p).hom

@[simp] theorem evalAffine_zero (A : CommRingCat.{0})
    (p : Spec (CommRingCat.of ℂ) ⟶ Spec A) : evalAffine A p 0 = 0 :=
  map_zero (evalHom A p).hom

@[simp] theorem evalAffine_add (A : CommRingCat.{0})
    (p : Spec (CommRingCat.of ℂ) ⟶ Spec A) (a b : A) :
    evalAffine A p (a + b) = evalAffine A p a + evalAffine A p b :=
  map_add (evalHom A p).hom a b

@[simp] theorem evalAffine_mul (A : CommRingCat.{0})
    (p : Spec (CommRingCat.of ℂ) ⟶ Spec A) (a b : A) :
    evalAffine A p (a * b) = evalAffine A p a * evalAffine A p b :=
  map_mul (evalHom A p).hom a b

/-! ## ★★★評価は点を決める(単射性) -/

/-- ★**`evalHom` は `Spec.map` の逆である**。 -/
@[simp] theorem Spec_map_evalHom (A : CommRingCat.{0})
    (p : Spec (CommRingCat.of ℂ) ⟶ Spec A) : Spec.map (evalHom A p) = p :=
  Spec.map_preimage p

/-- ★★★**評価は点を完全に決める**。

★★これは `topology_affine` が**意味を持つ**ための要である——
評価写像が点を潰していたら、そこから誘導した位相は点を分離できない。 -/
theorem evalHom_injective (A : CommRingCat.{0}) :
    Function.Injective (evalHom A) := by
  intro p q h
  rw [← Spec_map_evalHom A p, ← Spec_map_evalHom A q, h]

/-- ★**評価の値がすべて一致すれば同じ点である**。 -/
theorem eq_of_evalAffine_eq (A : CommRingCat.{0})
    {p q : Spec (CommRingCat.of ℂ) ⟶ Spec A}
    (h : ∀ a : A, evalAffine A p a = evalAffine A q a) : p = q := by
  refine evalHom_injective A ?_
  ext a
  exact h a

/-- ★★**逆に、環準同型からは点が作れる**——対応は全単射である。 -/
@[simp] theorem evalHom_Spec_map (A : CommRingCat.{0}) (φ : A ⟶ CommRingCat.of ℂ) :
    evalHom A (Spec.map φ) = φ :=
  Spec.preimage_map φ

/-! ## ★出典の紐付け(`.src`)

★条つきである。原文の `X^arc` は位相空間であり、
本ファイルはその**位相を固定する評価写像**だけを与える。 -/

def evalAffine.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——アフィンでの複素点の評価のみ)",
    sectionId := "genell-def-1-1-i" }

def evalHom_injective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——評価が点を決めること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
