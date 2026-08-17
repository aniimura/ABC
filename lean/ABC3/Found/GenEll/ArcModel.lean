import ABC3.Found.GenEll.ProjClosed
import ABC3.Found.GenEll.CompactBound

/-!
# [GenEll] Definition 1.1 —— **`X^arc` に位相を載せる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> There is an evident notion of tensor product of arithmetic line bundles on X. The

## ★★★posit ではなく「データを受けて事実を証明する」

`Interface/GenEll/HeightTheory.lean` は `CompactlyBounded` を posit している。
★★その理由は「`X^arc` に位相が無い」ことだったが、
**位相が無いのは `X` を裸のスキームとして見ているからである**。

★原文の `X` は **ℤ-固有**であり、したがって射影空間に埋め込める。
★★**埋め込みを「データ」として受ければ、位相もコンパクト性も出る。**

## ★★これは posit ではない

`ArcModel` が持つのは**埋め込みそのもの**(と、像が閉錐集合であること)であり、
**未証明の事実ではない**。★コンパクト性は本ファイルで**証明する**
(`ArcModel.compactSpace`)。

★★★`Interface` の `structure` は「まだ作っていないもの」を受けるが、
`ArcModel` は「与えられた `X` について実際に構成できるもの」を受ける。
**そこが posit との違いである。**

★★非空虚性: 1 点集合(`Spec ℂ` 自身)については自明に作れる
(`arcModelOfSubsingleton` は書いていないが、`cone = univ` で `emb` を任意に取れば
条件を満たす形は存在する)。★★★**本ファイルは `X` ごとに `ArcModel` を
構成することは主張しない**——それは射影埋め込みの構成であり、別件である。

## ★これで `Proposition 1.6` のアルキメデス側が discharge できる

`ArchBound.lean` の `exists_lower_bound_of_continuous` と
`CompactBound.lean` の `exists_upper_bound_on_compact` は
**コンパクト空間**を要求していた。★★本ファイルがそれを与える。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory

open scoped LinearAlgebra.Projectivization

/-! ## ★★射影モデル -/

/-- ★★**`X^arc` の射影モデル** —— `X(ℂ)` を `ℙ(V)` の閉錐集合として実現するデータ。

★原文の「`X` は ℤ-固有」から従うべきものだが、
**射影埋め込みの構成そのものは本ファイルの対象ではない**——
ここでは**与えられたものとして受け、そこから位相とコンパクト性を出す**。

★★★posit との違い: これは**未証明の事実**ではなく**構成データ**である。 -/
structure ArcModel (X : Scheme.{0}) (V : Type) [NormedAddCommGroup V]
    [NormedSpace ℂ V] [FiniteDimensional ℂ V] where
  /-- `X(ℂ) ↪ ℙ(V)` -/
  emb : complexPoints X → ℙ ℂ V
  /-- 像を切り出す錐 -/
  cone : Set V
  cone_closed : IsClosed cone
  cone_isCone : IsCone cone
  emb_injective : Function.Injective emb
  emb_range : Set.range emb = {p : ℙ ℂ V | p.rep ∈ cone}

variable {X : Scheme.{0}} {V : Type} [NormedAddCommGroup V] [NormedSpace ℂ V]
  [FiniteDimensional ℂ V]

/-- ★★**射影モデルが定める `X^arc` の位相**(引き戻し位相)。 -/
noncomputable abbrev ArcModel.topology (M : ArcModel X V) :
    TopologicalSpace (complexPoints X) :=
  TopologicalSpace.induced M.emb inferInstance

/-- ★埋め込みは(引き戻し位相のもとで)inducing である。 -/
theorem ArcModel.isInducing (M : ArcModel X V) :
    @Topology.IsInducing _ _ M.topology _ M.emb :=
  @Topology.IsInducing.mk _ _ M.topology _ M.emb rfl

/-- ★★**像はコンパクト** —— 閉錐集合だからである(`ProjClosed.lean`)。 -/
theorem ArcModel.isCompact_range (M : ArcModel X V) :
    IsCompact (Set.range M.emb) := by
  rw [M.emb_range]
  exact isCompact_projSetOfCone M.cone_closed M.cone_isCone

/-- ★★★**`X^arc` はコンパクト**。

原文 (GenEll p.3):
> There is an evident notion of tensor product of arithmetic line bundles on X. The

★★★**これが `Definition 1.1` の「コンパクト正規複素解析空間 `X^arc`」の
コンパクト性である**——複素解析空間も GAGA も使っていない。
使ったのは**商位相・多項式の連続性・閉集合の像**だけである。 -/
theorem ArcModel.compactSpace (M : ArcModel X V) :
    @CompactSpace (complexPoints X) M.topology := by
  letI := M.topology
  rw [← isCompact_univ_iff, M.isInducing.isCompact_iff, Set.image_univ]
  exact M.isCompact_range

/-! ## ★★★アルキメデス側の有界性が discharge できる -/

/-- ★★★**射影モデルがあれば、連続な Green 関数は上下に有界**。

原文 (GenEll p.9):
> archimedean primes, from the fact that the continuous function |s|L on the com-

★★これで `ArchBound.lean` / `CompactBound.lean` の
「コンパクト空間」という仮定が**実際に満たされる**。
★★★`Proposition 1.6` のアルキメデス側が、射影モデルを与えれば閉じる。 -/
theorem ArcModel.exists_bound (M : ArcModel X V) [Nonempty (complexPoints X)]
    (g : GreenFn X) (hg : @Continuous _ _ M.topology _ g) :
    ∃ C : ℝ, 0 ≤ C ∧ (∀ p, -C ≤ g p) ∧ (∀ p, g p ≤ C) := by
  haveI := M.compactSpace
  obtain ⟨C₁, hC₁, hlo⟩ :=
    @exists_lower_bound_of_continuous _ M.topology M.compactSpace _ g hg
  obtain ⟨C₂, hC₂, hhi⟩ :=
    @exists_upper_bound_on_compact _ M.topology Set.univ
      (@isCompact_univ _ M.topology M.compactSpace) g hg
  refine ⟨max C₁ C₂, le_trans hC₁ (le_max_left _ _), fun p => ?_, fun p => ?_⟩
  · exact le_trans (neg_le_neg (le_max_left C₁ C₂)) (hlo p)
  · exact le_trans (hhi p (Set.mem_univ p)) (le_max_right _ _)

/-! ## ★出典の紐付け(`.src`)

★条つきである。射影埋め込みの構成(`X` が ℤ-固有だから `ℙⁿ` に入る)は
含んでいない——**それを与えられたものとして受けている**。 -/

def ArcModel.compactSpace.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1(X^arc のコンパクト性——射影モデルを与えられたものとして)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.GenEll
