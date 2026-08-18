import ABC3.Found.Arakelov.PicGenSection

/-!
# Arakelov (B1) 第 114 ブロック —— **生成族だけで足りる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★「基底 vs すべて」の問題が消えた

第 101・102 ブロックは**被覆篩のすべての射**で全単射を要求していた。
★★しかし篩は**下方閉**なので、`D(t·g)` を入れるとその下の**すべての開集合**が入ってしまう
——そこでは全単射が言えない。

★★★**解決**: `imageSieve` も `equalizerSieve` も**篩**である。
したがって `Sieve.generate_le_iff` により、
**生成族(presieve)の射だけ**で全単射を言えば十分である。

## ★★これが §9-105 で特定した構造上の 1 点の**本当の解**である

★第 102 ブロックは器具としては正しかったが、仮定が**強すぎた**。
★★本ブロックで**生成族版**に弱めることで、基本開集合だけを見れば済むようになる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isLocallySurjective_of_presieve` | ★★★生成族で持ち上がれば局所全射 |
| `isLocallyInjective_of_presieve` | ★★★生成族で一致すれば局所単射 |
-/

universe v' u' v u

namespace ABC3.Found.Arakelov

open CategoryTheory Opposite

variable {C : Type u} [Category.{v} C] (J : GrothendieckTopology C)
  {A : Type u'} [Category.{v'} A] {FA : A → A → Type*} {CA : A → Type v'}
  [∀ X Y, FunLike (FA X Y) (CA X) (CA Y)] [ConcreteCategory.{v'} A FA]

/-- ★★★**生成族の上で持ち上がれば局所全射である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★`imageSieve` は篩なので、`Sieve.generate_le_iff` で生成族だけ見れば足りる。 -/
theorem isLocallySurjective_of_presieve {F G : Cᵒᵖ ⥤ A} (f : F ⟶ G)
    (h : ∀ (U : C) (s : ToType (G.obj (op U))), ∃ R : Presieve U, Sieve.generate R ∈ J U ∧
      ∀ ⦃W : C⦄ (i : W ⟶ U), R i → ∃ t, f.app (op W) t = G.map i.op s) :
    Presheaf.IsLocallySurjective J f where
  imageSieve_mem {U} s := by
    obtain ⟨R, hR, hb⟩ := h U s
    exact J.superset_covering ((Sieve.generate_le_iff R _).2 (fun {W} i hi => hb i hi)) hR

/-- ★★★**生成族の上で一致すれば局所単射である**。 -/
theorem isLocallyInjective_of_presieve {F G : Cᵒᵖ ⥤ A} (f : F ⟶ G)
    (h : ∀ (U : C) (x y : ToType (F.obj (op U))), f.app (op U) x = f.app (op U) y →
      ∃ R : Presieve U, Sieve.generate R ∈ J U ∧
        ∀ ⦃W : C⦄ (i : W ⟶ U), R i → F.map i.op x = F.map i.op y) :
    Presheaf.IsLocallyInjective J f where
  equalizerSieve_mem {X} x y hxy := by
    obtain ⟨R, hR, hb⟩ := h X.unop x y hxy
    exact J.superset_covering ((Sieve.generate_le_iff R _).2 (fun {W} i hi => hb i hi)) hR

/-! ## ★出典の紐付け(`.src`) -/

def isLocallySurjective_of_presieve.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——生成族だけで局所全単射が出る)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
