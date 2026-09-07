import ABC3.Meta.Claim
import Mathlib.FieldTheory.KrullTopology
import Mathlib.NumberTheory.Cyclotomic.CyclotomicCharacter
import Mathlib.NumberTheory.Padics.PadicNumbers
import Mathlib.FieldTheory.IsAlgClosed.AlgebraicClosure
import Mathlib.Topology.Algebra.ContinuousMonoidHom

/-!
# [pGC] §1 の設定 — 記号と、地の文に埋め込まれた定義

原典: S. Mochizuki, *A Version of the Grothendieck Conjecture for p-adic Local Fields* (1997)。
構造化: `ResearchPaper/1_Structured/A Version of the Grothendieck Conjecture for p-adic Local Fields/section-1.html`
(PDF 目視確認済み: 物理 p.2・p.3・p.4冒頭)。

## この段階で確定していること / していないこと

- 型は付く。**証明は付けていない**(`sorry`)。
- §1 冒頭の「group-theoretically に回復できる」は、原典が**見出し語を持たない地の文**で
  与えている定義。ここで明示的に切り出した(`RecoverableFromAbsGal`)。
- ★**原典が暗黙にしていて、我々が明示にしたもの**: 対象を α に沿って移す規則(`transport`)。
  原典は「the object associated to K is necessarily taken by α to the corresponding object
  associated to K′」と述べるだけで、移送の規則自体は文脈から自明として与えていない。
  これをデータとして持たせたのは**我々の設計判断**であり、原典の転写ではない。
  移送の選び方によって主張の強さが変わるため、下記の識別力検査を併置する。
-/

namespace ABC3.Skeleton.PGC

open ABC3.Meta

variable (p : ℕ) [Fact p.Prime]

/-! ## 設定 -/

/-- 本論文でいう **p進局所体** — `ℚ_[p]` の有限次拡大。

原典 [pGC] 物理 p.2:
「Let p be a prime number. Let K be a p-adic local field. By this, we shall mean in this
paper that K is a finite extension of **Q**_p.」 -/
structure PAdicLocalField where
  carrier : Type
  [isField : Field carrier]
  [isAlgebra : Algebra ℚ_[p] carrier]
  [isFinite : FiniteDimensional ℚ_[p] carrier]

def PAdicLocalField.src : Source :=
  { paper := "pGC", pdfPage := 2, item := "Section 1 (opening)", sectionId := "setup-1-intro" }

attribute [instance] PAdicLocalField.isField PAdicLocalField.isAlgebra PAdicLocalField.isFinite

namespace PAdicLocalField

variable {p}

/-- 原典の `K̄` — K の代数閉包(原典は「Fix an algebraic closure K̄ of K」)。 -/
abbrev closure (K : PAdicLocalField p) : Type := AlgebraicClosure K.carrier

def closure.src : Source :=
  { paper := "pGC", pdfPage := 2, item := "Section 1 (opening)", sectionId := "setup-1-intro" }

/-- 原典の `Γ_K ≝ Gal(K̄/K)` — 絶対 Galois 群。

位相は mathlib の `krullTopology` が与える(`Mathlib/FieldTheory/KrullTopology.lean`)。
原典は「isomorphism of topological groups」と言うので、位相込みで扱う必要がある。 -/
abbrev absGal (K : PAdicLocalField p) : Type := K.closure ≃ₐ[K.carrier] K.closure

def absGal.src : Source :=
  { paper := "pGC", pdfPage := 2, item := "Section 1 (opening)", sectionId := "setup-1-intro" }

end PAdicLocalField

/-! ## §1 冒頭の暗黙の定義 -/

/-- **K に付随する対象**の族と、Γ_K の同型に沿った移送の規則。

原典は「an object associated to K」という言い方をするだけで、対象の種類を限定していない
——実際 §1 では3種類が現れる: 指標 χ(Γ_K 上の関数)・数 q や [K:ℚ_p](Γ_K に依存しない値)・
惰性群 I_K(Γ_K の部分群)。ゆえに対象の**型そのもの**が K ごとに変わりうる形にしてある。

★`transport` は原典が明示していない(上記モジュール docstring 参照)。 -/
structure AssociatedObject where
  /-- K に付随する対象の型 -/
  Obj : PAdicLocalField p → Type
  /-- K に付随する対象そのもの -/
  obj : (K : PAdicLocalField p) → Obj K
  /-- 位相群の同型 α : Γ_K ≅ Γ_K′ に沿った移送 -/
  transport : {K K' : PAdicLocalField p} →
    ContinuousMulEquiv K.absGal K'.absGal → Obj K → Obj K'

def AssociatedObject.src : Source :=
  { paper := "pGC", pdfPage := 2, item := "Section 1 (opening)", sectionId := "setup-1-intro" }

variable {p}

/-- **「Γ_K から group-theoretically に回復できる」** — [pGC] §1 冒頭の暗黙の定義。

原文 (pGC p.2):
> when we say that an object associated to K can be recovered "group-theoretically" from
> Γ_K, we mean that given another local p-adic field K[prime], together with an isomorphism
> of topological groups α : Γ_K ≅ Γ[prime]_K, the object associated to K is necessarily
> taken by α to the corresponding object associated to K[prime].

注意すべき形: これは**単称の構成手続きの存在**ではなく、**すべての同型に対する両立性**である。
「Γ_K から χ を作るアルゴリズムがある」ではなく「どの同型で写しても対応が保たれる」。 -/
def AssociatedObject.RecoverableFromAbsGal (A : AssociatedObject p) : Prop :=
  ∀ {K K' : PAdicLocalField p} (α : ContinuousMulEquiv K.absGal K'.absGal),
    A.transport α (A.obj K) = A.obj K'

def AssociatedObject.RecoverableFromAbsGal.src : Source :=
  { paper := "pGC", pdfPage := 2, item := "Section 1 (opening)", sectionId := "setup-1-intro" }

/-!  ## 識別力の検査

`RecoverableFromAbsGal` が常に真でも常に偽でもないことの検査は、主語が違う
(原典の主張ではなく**我々のモデル**についての事実)ため `ABC3/Check/PGC/Section1Discriminating.lean`
に置いた。`Skeleton/` から `Check/` を import してはならない。
-/

/-! ## [pGC] Definition 2.3 — filtered group

原典 [pGC] 物理 p.5:
> We shall call a collection of closed normal subgroups {Gᵛ} of G (where v ranges over all
> positive real numbers) a filtration on G if Gᵛ¹ ⊆ Gᵛ² whenever v1 ≥ v2.

構造化: `ResearchPaper/1_Structured/A Version of the Grothendieck Conjecture for p-adic Local
Fields/section-2.html#def-2-3`(PDF 目視確認 2026-09-03、物理 p.5)。

## ★なぜ `Found/PGC/FilteredGroup.lean` からここへ移したか(2026-09-08)

移動前、この定義は `Found/PGC/FilteredGroup.lean` に**だけ**あった。ところが
`Skeleton/PGC/Section3.lean`(`filteredGroupOf`)と `Skeleton/PGC/Section4.lean`
(`RamificationFiltration.filt`)は**定義と定理が同居した本**であり、その定義を
`Found` から引くには `Skeleton/PGC/Section3` ごと import するしかなかった。
`node tools/import-audit.mjs`(2026-09-07 実測)によると、そのせいで
`Found/PGC/` の 202 本のうち **`Section3` は 200 本・`Section4` は 202 本**が
「import すると循環する」状態だった。★詰まる原因は 1 つ、同じもの——
**`FilteredGroup` が `Found` の本にしかない**こと。

`Section1Defs.lean`・`Section2Defs.lean` と同じ作法で定義を割ろうとしても、
割った先の `Section3Defs` / `Section4Defs` が結局
`ABC3.Found.PGC.FilteredGroup` を持ち越す——つまり **Skeleton の Defs が
Found を import する**ことになり、割る作業が 1 段で閉じない。
そこで定義そのものをここ(`Skeleton`)へ降ろした。

## ★なぜ `Interface/PGC/` ではなくここか(測って決めた)

`Interface/PGC/LocalFieldData.lean` には `ResidueCardinality` /
`RamificationFiltration` という同種の `structure` があるので候補ではあったが、
次の 2 点を測って `Skeleton/PGC/Setup.lean` を選んだ:

1. `Interface/` の `structure` は `tools/check.mjs` の **G2**(実装は
   `check.mjs:1041`)により `X.nonvacuous` か `X.waiting` を要求される。
   `FilteredGroup` と `FilteredGroup.Iso` は `structure` なので **witness が
   2 つ増える**。しかも `LocalFieldData.lean` 冒頭の規約により witness は
   `Found/` 側に置かねばならず、新しい本が 1 つ増える。
2. `Interface/` の定義は「原典が所与として使うが、我々がまだ mathlib から
   得られないもの」である(同ファイル冒頭)。`FilteredGroup` は
   `sorry` 無しで完全に構成済みで、境界外入力を 1 つも要求しない
   ——`Interface` の役割に当てはまらない。

一方 `Skeleton/PGC/Setup.lean` は「記号と、地の文に埋め込まれた定義」の置き場で、
`Interface` / `Skeleton` / `Found` の PGC 側**すべて**が既に(推移的に)
import している。ゆえに**新しい import 辺は 1 本も増えない**。

★`structure` は `check.mjs` の G8/G9 の対象外(両ゲートは `theorem`/`lemma` のみを見る)
なので、「実装は `Found/` に置く」という規約には触れない。

★旧名 `ABC3.Found.PGC.FilteredGroup` は `Found/PGC/FilteredGroup.lean` に残した
`export` で引き続き解決できる(そちらの docstring を参照)。

## 逸脱の記録(CLAUDE.md「逸脱」)

原文はこの直後で「射・同型・外部同型の定義は読者に委ねる」と明言しており(§2 の
「省略」の合図そのもの)、下記の `FilteredGroup.Iso` / `FilteredGroup.OuterIso` は
**我々自身の定式化**である:

* **射**: 添字を保つ連続群射(`f (Gᵛ) ⊆ G'ᵛ`、同じ v で比較する)。
  原文は添字の対応規則を指定していないので、最も忠実な最小限の選択として
  「同じ添字」を採用した。
* **同型**: 添字を保つ連続同型で、`Gᵛ` の像がちょうど `G'ᵛ` に一致するもの
  (単なる `⊆` ではなく `=`——可逆なので両方向の `⊆` と同値)。
* **外部同型**: 同型を、後合成による内部自己同型の作用で割った商。
-/

/-- **[pGC] Definition 2.3** — filtered group。

`G` 上の実数で添字付けられた閉正規部分群の族 `Gv : ℝ → Subgroup G` であって、
添字が大きいほど小さい部分群になっている(降下条件)もの。

原文の条件「`v1 ≥ v2` ⟹ `Gᵛ¹ ⊆ Gᵛ²`」は `Antitone Gv`(`≤` の通常の順序に関する
反単調性)そのもの——`Antitone f` の定義 `a ≤ b → f b ≤ f a` に `a := v2, b := v1` を
代入すれば原文の条件と一致する。 -/
structure FilteredGroup where
  /-- 台となる位相群 -/
  G : Type*
  [isGroup : Group G]
  [isTop : TopologicalSpace G]
  [isTopGroup : IsTopologicalGroup G]
  /-- フィルトレーション `{Gᵛ}_{v>0}`。原文は `v` の範囲を正の実数に限るが、
  `v ≤ 0` での値を固定しても一般性は失わない(`Antitone` の条件は全実数上で書ける)。 -/
  Gv : ℝ → Subgroup G
  /-- 各 `Gᵛ` は閉部分群 -/
  isClosed : ∀ v, IsClosed (Gv v : Set G)
  /-- 各 `Gᵛ` は正規部分群 -/
  isNormal : ∀ v, (Gv v).Normal
  /-- 降下条件: `v1 ≥ v2 → Gᵛ¹ ⊆ Gᵛ²` -/
  antitone : Antitone Gv

attribute [instance] FilteredGroup.isGroup FilteredGroup.isTop FilteredGroup.isTopGroup

def FilteredGroup.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Definition 2.3", sectionId := "def-2-3" }

/-- filtered group の同型: 添字を保つ連続同型。 -/
structure FilteredGroup.Iso (A B : FilteredGroup) where
  /-- 台となる連続群同型 -/
  equiv : ContinuousMulEquiv A.G B.G
  /-- 各 v でフィルトレーションを保つ(像がちょうど一致する) -/
  map_Gv : ∀ v, Subgroup.map equiv.toMulEquiv (A.Gv v) = B.Gv v

/-- 原文は「射・同型の定義そのものは読者に委ねる」と明言しているので、
これは**我々自身の定式化**である(bare な `"Definition 2.3"` と紛れないよう、
item 名に注記を付ける——`RecoverableAsAddModule.src` と同じ作法)。 -/
def FilteredGroup.Iso.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Definition 2.3 (FilteredGroup.Iso)",
    sectionId := "def-2-3" }

/-- 内部自己同型で割った同値関係: `f ~ g` ⟺ ある `c : B.G` で `g = (c·(-)·c⁻¹) ∘ f`。 -/
instance FilteredGroup.Iso.setoid (A B : FilteredGroup) : Setoid (FilteredGroup.Iso A B) where
  r f g := ∃ c : B.G, ∀ x, g.equiv x = c * f.equiv x * c⁻¹
  iseqv := by
    refine ⟨fun f => ⟨1, by simp⟩, ?_, ?_⟩
    · rintro f g ⟨c, hc⟩
      exact ⟨c⁻¹, fun x => by rw [hc x]; group⟩
    · rintro f g h ⟨c, hc⟩ ⟨d, hd⟩
      exact ⟨d * c, fun x => by rw [hd, hc]; group⟩

/-- **外部同型**の集合: `FilteredGroup.Iso A B` を内部自己同型で割った商。

原文 [pGC] 物理 p.7 の `OutFilt(Γ_K, Γ_K')`(Theorem 4.2)はこれ。 -/
def FilteredGroup.OuterIso (A B : FilteredGroup) : Type _ :=
  Quotient (inferInstance : Setoid (FilteredGroup.Iso A B))

/-- 原文 p.5 は外部同型の定義も読者に委ねている——これも**我々自身の定式化**
(原文 p.7 の `OutFilt(Γ_K, Γ_K')` がこれに当たる)。 -/
def FilteredGroup.OuterIso.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Definition 2.3 (FilteredGroup.OuterIso)",
    sectionId := "def-2-3" }

end ABC3.Skeleton.PGC
