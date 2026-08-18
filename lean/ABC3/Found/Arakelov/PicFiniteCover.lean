import ABC3.Found.Arakelov.PicBasicTrivial

/-!
# Arakelov (B1) 第 135 ブロック —— **有限被覆と単位イデアル**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★`f^n` 論法には**有限性**が要る

古典的な筋(Hartshorne II.5.1)の 2 つの段

    (a) s|_{D f} = 0 ⟹ ∃ n, fⁿ s = 0
    (b) t ∈ F(D f)  ⟹ ∃ n, fⁿ t が ⊤ へ延びる

はどちらも「各 `D(gᵢ)` で得た `nᵢ` の**最大値**を取る」ので、
★★被覆が**有限**でなければならない。

## ★★取り方は「1 が有限個で書ける」

`Spec R` の準コンパクト性を経由せずとも、

    ⨆ D(G x) = ⊤  ⟺  Ideal.span (range G) = ⊤  ⟹  1 ∈ span (有限部分集合)

で有限部分被覆が出る(`Submodule.mem_span_finite_of_mem_span`)。
★添字を `Fin n` にしておくと後の `Finset.sup` が楽である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_finite_basicOpen_trivial` | ★★★**有限個の基本開集合で覆い、各々で自明** |

## ★★★新しい逃げ道——`beta_reduce`

`refine ⟨n, fun i => …, ?_, ?_⟩` の後、ゴールに `(fun i => …) i` が残って
`rw`/`▸` が当たらない。★`beta_reduce` 一発で通る。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

/-- ★★★**有限個の基本開集合で覆い、各々で自明化できる**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★`Ideal.span (range g) = ⊤` が「被覆」の代数的な言い換えである。 -/
theorem exists_finite_basicOpen_trivial (R : CommRingCat.{u}) (P : (Spec R).PresheafOfModules)
    (h : IsLocallyTrivial (Spec R) P) :
    ∃ (n : ℕ) (g : Fin n → (R : Type u)),
      Ideal.span (Set.range g) = ⊤ ∧
      ∀ i, Nonempty ((restrictPresheafFunctor (Spec R) (PrimeSpectrum.basicOpen (g i))).obj P
        ≅ 𝟙_ (PresheafModulesOn (Spec R) (PrimeSpectrum.basicOpen (g i)))) := by
  choose G hGmem hGtriv using exists_basicOpen_trivial R P h
  have hsup : (⨆ x : (Spec R), PrimeSpectrum.basicOpen (G x)) = ⊤ := by
    refine eq_top_iff.2 fun x _ => ?_
    exact Opens.mem_iSup.2 ⟨x, hGmem x⟩
  have hspan : Ideal.span (Set.range G) = ⊤ :=
    PrimeSpectrum.iSup_basicOpen_eq_top_iff.1 hsup
  obtain ⟨T, hTsub, hT1⟩ := Submodule.mem_span_finite_of_mem_span
    (show (1 : (R : Type u)) ∈ Ideal.span (Set.range G) by rw [hspan]; trivial)
  refine ⟨T.card, fun i => (T.equivFin.symm i : (R : Type u)), ?_, fun i => ?_⟩
  · refine Ideal.eq_top_iff_one _ |>.2 ?_
    have : Set.range (fun i => (T.equivFin.symm i : (R : Type u))) = (T : Set (R : Type u)) := by
      ext y
      constructor
      · rintro ⟨i, rfl⟩; exact (T.equivFin.symm i).2
      · intro hy; exact ⟨T.equivFin ⟨y, hy⟩, by simp⟩
    rw [this]
    exact hT1
  · obtain ⟨x, hx⟩ := hTsub (T.equivFin.symm i).2
    beta_reduce
    exact hx ▸ hGtriv x


/-! ## ★出典の紐付け(`.src`) -/

def exists_finite_basicOpen_trivial.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——有限被覆と単位イデアル)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
