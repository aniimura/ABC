import ABC3.Found.Arakelov.PicAwayFreeIso

/-!
# Arakelov (B1) 第 101 ブロック —— **被覆から局所全単射を出す器具**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★`Over V` の site で「基底 ⟹ 全体」を言うための土台

§9-105 で特定した構造上の 1 点——`IsLocallyTrivial` は
`Over V` **上の前層加群としての同型**を要求するが、
切断ごとの同型は**基底の上でしか**出ていない。

★★§9-106 で方針を決めた: `Sheaf.isLocallyBijective_iff_isIso`
(局所全単射 ⟺ 同型)は **`over` 位相でも使える**(2026-08-18 実測——
§9-73 で詰まったのは `SheafOfModules.toSheaf` 固有の問題だった)。

★★★あとは「基底で同型」から**局所全単射**を出せば良い。本ブロックがその器具である。

## ★★機構 —— 被覆篩を `imageSieve` / `equalizerSieve` に含める

★`imageSieve f s` は「`s` の制限が像に入る射」の篩、
`equalizerSieve x y` は「制限すると一致する射」の篩である。
★★被覆篩がそれらに**含まれれば**、`J.superset_covering` で被覆になる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isLocallySurjective_of_cover` | ★★★**被覆の上で持ち上がれば局所全射** |
| `isLocallyInjective_of_coverSieve` | ★★★**被覆の上で一致すれば局所単射** |
-/

universe v' u' v u

namespace ABC3.Found.Arakelov

open CategoryTheory Opposite

variable {C : Type u} [Category.{v} C] (J : GrothendieckTopology C)
  {A : Type u'} [Category.{v'} A] {FA : A → A → Type*} {CA : A → Type v'}
  [∀ X Y, FunLike (FA X Y) (CA X) (CA Y)] [ConcreteCategory.{v'} A FA]

/-- ★★★**被覆の上で持ち上がれば局所全射である**。 -/
theorem isLocallySurjective_of_cover {F G : Cᵒᵖ ⥤ A} (f : F ⟶ G)
    (h : ∀ (U : C) (s : ToType (G.obj (op U))), ∃ S : Sieve U, S ∈ J U ∧
      ∀ ⦃W : C⦄ (i : W ⟶ U), S.arrows i → ∃ t, f.app (op W) t = G.map i.op s) :
    Presheaf.IsLocallySurjective J f where
  imageSieve_mem {U} s := by
    obtain ⟨S, hS, hsub⟩ := h U s
    exact J.superset_covering (fun {W} i hi => hsub i hi) hS

/-- ★★★**被覆の上で一致すれば局所単射である**。 -/
theorem isLocallyInjective_of_coverSieve {F G : Cᵒᵖ ⥤ A} (f : F ⟶ G)
    (h : ∀ (U : C) (x y : ToType (F.obj (op U))), f.app (op U) x = f.app (op U) y →
      ∃ S : Sieve U, S ∈ J U ∧
        ∀ ⦃W : C⦄ (i : W ⟶ U), S.arrows i → F.map i.op x = F.map i.op y) :
    Presheaf.IsLocallyInjective J f where
  equalizerSieve_mem {X} x y hxy := by
    obtain ⟨S, hS, hsub⟩ := h X.unop x y hxy
    exact J.superset_covering (fun {W} i hi => hsub i hi) hS

/-! ## ★出典の紐付け(`.src`) -/

def isLocallySurjective_of_cover.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——被覆から局所全単射を出す器具)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
