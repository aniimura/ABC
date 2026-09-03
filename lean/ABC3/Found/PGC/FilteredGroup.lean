import ABC3.Skeleton.PGC.Setup

/-!
# [pGC] Definition 2.3 — filtered group

原典 [pGC] 物理 p.5:
> We shall call a collection of closed normal subgroups {Gᵛ} of G (where v ranges over all
> positive real numbers) a filtration on G if Gᵛ¹ ⊆ Gᵛ² whenever v1 ≥ v2.

構造化: `ResearchPaper/1_Structured/A Version of the Grothendieck Conjecture for p-adic Local
Fields/section-2.html#def-2-3`(PDF 目視確認 2026-09-03、物理 p.5)。

## なぜ `Found/` に直接置くか

Definition 2.3 は**純粋な定義**であり、原典の他の項目(局所 Tate 双対性・局所類体論・
高次分岐群の Herbrand の定理・p 進対数)のような境界外入力を一切要求しない——`G` は
**任意の**位相群でよく、`Γ_K` に固有の性質は使わない。ゆえに `Skeleton/` を経由せず、
**`sorry` 無しでここに直接置く**(`structure` は `tools/check.mjs` の G8/G9 の対象外
——両ゲートは `theorem`/`lemma` のみを見る)。

## 逸脱の記録(CLAUDE.md「逸脱」)

原文はこの直後で「射・同型・外部同型の定義は読者に委ねる」と明言しており(§2 の
「省略」の合図そのもの)、これらは本ファイルでは扱わない——`Skeleton/PGC/Section2.lean`
の docstring に、その3つの定義の設計判断を記録する。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC

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

/-! ## filtered group の射・同型・外部同型

原文 [pGC] 物理 p.5 は「射・同型の定義そのものは読者に委ねる」と明言している
(`Skeleton/PGC/Section2.lean` の `remark-2-3-morphisms` 読みを参照)——本節がここで
定式化するのは我々自身の判断であり、**逸脱として記録する**:

* **射**: 添字を保つ連続群射(`f (Gᵛ) ⊆ G'ᵛ`、同じ v で比較する)。
  原文は添字の対応規則を指定していないので、最も忠実な最小限の選択として
  「同じ添字」を採用した。
* **同型**: 添字を保つ連続同型で、`Gᵛ` の像がちょうど `G'ᵛ` に一致するもの
  (単なる `⊆` ではなく `=`——可逆なので両方向の `⊆` と同値)。
* **外部同型**: 同型を、後合成による内部自己同型の作用で割った商。 -/

/-- filtered group の同型: 添字を保つ連続同型。 -/
structure FilteredGroup.Iso (A B : FilteredGroup) where
  /-- 台となる連続群同型 -/
  equiv : ContinuousMulEquiv A.G B.G
  /-- 各 v でフィルトレーションを保つ(像がちょうど一致する) -/
  map_Gv : ∀ v, Subgroup.map equiv.toMulEquiv (A.Gv v) = B.Gv v

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

end ABC3.Found.PGC
