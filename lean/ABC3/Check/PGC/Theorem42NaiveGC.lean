import ABC3.Check.PGC.Prop12ForallRD
import ABC3.Found.PGC.RamificationNaturality

/-!
# [pGC] Theorem 4.2 の**現行形**は「素朴な Grothendieck 予想」を含む

原文 (pGC p.7, Theorem 4.2):

> Let K and K[prime] be finite extensions of Q_p. Let Isom_Q_p(K, K[prime]) denote the set of
> Q_p-algebra isomorphisms of K with K[prime]. Let Out_Filt(Γ_K, Γ_K[prime]) denote the set of outer
> isomorphisms of filtered groups between the absolute Galois groups of K and K[prime] equipped
> with the filtrations defined by the higher (i.e., with index > 0) ramification groups in
> the upper numbering. Then the natural morphism
> Isom_Q_p(K, K[prime]) → Out_Filt(Γ_K, Γ_K[prime])
> induced by "looking at the morphism induced on absolute Galois groups" is a bijection.

`Skeleton/PGC/Section4.lean:122` の現行形は

```
theorem theorem_4_2 (RF : RamificationFiltration p) (hnat : IsNaturalFiltration RF)
    (K K' : PAdicLocalField p) :
    Function.Bijective (naturalOuterIso RF hnat (K := K) (K' := K')) := sorry
```

である。2026-09-05 に「`Φ` を自由な関数で受け取っていた旧形(`Check/PGC/Theorem42Degenerate.lean`
で反例つき)を `IsNaturalFiltration` を足して修理した」と記録されている。

## 修理は効いていない——`IsNaturalFiltration` は何も切らない

`Found/PGC/RamificationNaturality.lean:75` の `topFiltration`(`Gv _ _ := ⊤`)が
`IsNaturalFiltration` を満たすことは同ファイル `:81` の `exists_isNaturalFiltration` で
2 行で証明済みである(本ファイルではその証明本体を `topFiltration_isNatural` として
名前付きで取り出す)。ゆえに `RF := topFiltration p` を代入できる。

代入すると `FilteredGroup.Iso` の `map_Gv` 条件は `Subgroup.map φ ⊤ = ⊤` になり、
これは `φ` の全射性だけから従う——`φ` に**何の制約も課さない**。すなわち

```
FilteredGroup.OuterIso (filtOf (topFiltration p) K) (filtOf (topFiltration p) K')
  =  Γ_K と Γ_{K'} の連続同型の全体 / 内部自己同型
```

であり、これは原文が「filtration との両立を外した」外部同型の全体そのものである
(`nonempty_topOuterIso_iff` が非空性のレベルでこれを述べる)。

したがって現行形の**全射性**は「任意の連続外部同型が体の同型から来る」を主張する。
それが `theorem_4_2_current_form_implies_naive_GC`(`sorry` 無し)である。
★単射性は一切使わない(`surjective_naturalOuterIso_implies_naive_GC`)。

## なぜそれが問題なのか

右辺は原典が Introduction で明示的に偽だと述べている命題である。

原文 (pGC p.1, Introduction):

> On the one hand, one knows (cf. the Remark in [4] following Theorem 4.2) that the
> Grothendieck Conjecture cannot hold in the naive sense (i.e., if one removes the condition
> of "compatibility with the filtrations" from the outer isomorphisms considered - see, e.g.,
> [8]), so one must put some sort of condition on the outer isomorphisms of Galois groups that
> one considers.

原文 (pGC p.1, Historical Remark):

> I originally set out to prove the naive version of the above Theorem, only to discover
> that this was, in fact, false.

つまり現行形は、まさに原典が「これは偽だから条件を課さねばならない」と言っている
「条件を外した版」を(特殊化として)含んでいる。`IsNaturalFiltration` は
「射 `naturalOuterIso` を**構成する**ために必要な仮説」であって、
「フィルトレーションが本物の高次分岐群であること」を言う仮説ではない——
そこが `Interface.PGC.RamificationFiltration` の抽象性ごと抜けている。

## 「落とした条件は主張を偽にするか自明にする」——修理が穴を移しただけの例

旧形の穴は「射 `Φ` が自由」だった。現行形は射を構成して塞いだが、
**フィルトレーション `RF` の側が自由なまま**なので、`RF` を退化させれば
同じ穴から同じ結論(原典が偽と述べた命題)が出る。第 1013
(`Check/PGC/Prop12ForallRD.lean`)が Proposition 1.2 について見つけたのと
同じ穴であり、`theorem_4_2_current_form_implies_forall_RD` で実際に繋げてある。

## 正しい形の見通し(本ファイルでは実装しない)

`RF` を全称量化から外し、本物の高次分岐群(上付き番号付け)として構成された
唯一の `RF` に固定する必要がある。それには Herbrand の定理が要る
(`Skeleton/PGC/Section4.lean` の `.needs`、`memory/pgc-ramification-naturality-gap.md`)。
★`Skeleton` の statement をどう直すかは人の判断待ち(`decisions-pending.md` D24)。
本ファイルは判断の材料を出すだけで、`Skeleton`/`Interface` には触れていない。

**これは原典の主張ではない**(我々のモデルと器具についての事実)ので `.src` を持たない。

## 逸脱の記録(CLAUDE.md「逸脱」)

* 仮説として書いた `∀ RF hnat K K', Function.Bijective (naturalOuterIso RF hnat)` は
  `Skeleton/PGC/Section4.lean::theorem_4_2` の statement を全称量化しただけのものだが、
  `Skeleton` を import せずに済むよう**型を写して**いる(`theorem_4_2` 自体は `sorry`
  なので、それを使うと `sorryAx` が axiom 一覧に入ってしまう)。写しが原文と一致することは
  上のコードブロックで目視できる。
* `topFiltration_isNatural` は `Found/PGC/RamificationNaturality.lean` に置くのが
  筋だが、そうすると `Skeleton/PGC/Section4` 以降が全て再ビルドになるため、
  本ファイル内に置いた(証明は `exists_isNaturalFiltration` の本体と同一)。
-/

namespace ABC3.Check.PGC

open ABC3.Skeleton.PGC ABC3.Found.PGC ABC3.Interface.PGC

variable (p : ℕ) [Fact p.Prime]

/-! ## 1. `topFiltration` は `IsNaturalFiltration` を満たす -/

/-- 退化フィルトレーション `Gv ≡ ⊤` は `IsNaturalFiltration` を満たす。

`Found/PGC/RamificationNaturality.lean::exists_isNaturalFiltration` の証明本体を
そのまま名前付きの補題として取り出したもの(あちらは `∃` の形でしか公開していない)。
`RF := topFiltration p` を代入するために必要。 -/
theorem topFiltration_isNatural : IsNaturalFiltration (topFiltration p) := by
  intro K K' α v
  exact Subgroup.map_top_of_surjective _ (galContinuousMulEquiv α).surjective

/-! ## 2. `topFiltration` を入れると `OuterIso` は「連続外部同型の全体」になる -/

/-- 連続同型 `φ : Γ_K ≃ₜ* Γ_{K'}` から、退化フィルトレーション付きの
filtered group の同型を作る。

`map_Gv` の中身は `Subgroup.map φ ⊤ = ⊤`——すなわち `φ` の全射性だけで従い、
`φ` に何の制約も課さない。これが「`topFiltration` を入れると条件が消える」の中身。 -/
noncomputable def topIsoOfContinuousMulEquiv {K K' : PAdicLocalField p}
    (φ : ContinuousMulEquiv K.absGal K'.absGal) :
    FilteredGroup.Iso (filtOf (topFiltration p) K) (filtOf (topFiltration p) K') where
  equiv := φ
  map_Gv _ := Subgroup.map_top_of_surjective _ φ.surjective

/-- **★★★★★`topFiltration` 代入後の `OuterIso` が実際に何になったか**。

`Out_Filt(Γ_K, Γ_{K'})` は `Γ_K` と `Γ_{K'}` の連続同型の全体を内部自己同型で
割ったもの——すなわち原文が言う「filtration との両立条件を外した」外部同型の全体
そのものである。非空性のレベルでこれを述べる。 -/
theorem nonempty_topOuterIso_iff (K K' : PAdicLocalField p) :
    Nonempty (FilteredGroup.OuterIso (filtOf (topFiltration p) K) (filtOf (topFiltration p) K'))
      ↔ Nonempty (ContinuousMulEquiv K.absGal K'.absGal) := by
  constructor
  · rintro ⟨c⟩
    obtain ⟨f, -⟩ := Quotient.exists_rep c
    exact ⟨f.equiv⟩
  · rintro ⟨φ⟩
    exact ⟨Quotient.mk _ (topIsoOfContinuousMulEquiv p φ)⟩

/-! ## 3. 主定理——現行形は素朴 Grothendieck 予想を含む -/

/-- 使うのは**全射性だけ**である、を明示した形。

`RF := topFiltration p` を代入し、`φ` から作った外部同型類の原像を取るだけ。
単射性(`Function.Bijective` の第 1 成分)は一切使わない。 -/
theorem surjective_naturalOuterIso_implies_naive_GC
    (h : ∀ (RF : RamificationFiltration p) (hnat : IsNaturalFiltration RF)
      (K K' : PAdicLocalField p),
      Function.Surjective (naturalOuterIso RF hnat (K := K) (K' := K'))) :
    ∀ (K K' : PAdicLocalField p), ContinuousMulEquiv K.absGal K'.absGal →
      Nonempty (K.carrier ≃ₐ[ℚ_[p]] K'.carrier) := by
  intro K K' φ
  obtain ⟨α, -⟩ := h (topFiltration p) (topFiltration_isNatural p) K K'
    (Quotient.mk _ (topIsoOfContinuousMulEquiv p φ))
  exact ⟨α⟩

/-- **★★★★★★★★★★★★★★★★★★★★現行形の `theorem_4_2` は素朴 Grothendieck 予想を含む**。

仮説は `Skeleton/PGC/Section4.lean::theorem_4_2` の statement を
`RF`・`hnat`・`K`・`K'` について全称量化しただけのもの
(`Skeleton` を import せずに済むよう、型だけをここに写した——上の逸脱の記録を参照)。 -/
theorem theorem_4_2_current_form_implies_naive_GC
    (h : ∀ (RF : RamificationFiltration p) (hnat : IsNaturalFiltration RF)
      (K K' : PAdicLocalField p),
      Function.Bijective (naturalOuterIso RF hnat (K := K) (K' := K'))) :
    ∀ (K K' : PAdicLocalField p), ContinuousMulEquiv K.absGal K'.absGal →
      Nonempty (K.carrier ≃ₐ[ℚ_[p]] K'.carrier) :=
  surjective_naturalOuterIso_implies_naive_GC p
    (fun RF hnat K K' => (h RF hnat K K').2)

#print axioms theorem_4_2_current_form_implies_naive_GC

/-! ## 4. 帰結 -/

/-- 現行形は `Check/PGC/Prop12ForallRD.lean` が「原典が偽と述べている命題と同値」と
判定した `∀ RD` 版をも含む。すなわち D24 と第 1013 の欠陥は同じ穴である。 -/
theorem theorem_4_2_current_form_implies_forall_RD
    (h : ∀ (RF : RamificationFiltration p) (hnat : IsNaturalFiltration RF)
      (K K' : PAdicLocalField p),
      Function.Bijective (naturalOuterIso RF hnat (K := K) (K' := K'))) :
    ∀ RD : ResidueCardinality p, (residueCardAndDegreeObject RD).RecoverableFromAbsGal :=
  (forall_RD_recoverable_iff_algEquiv p).mpr
    (fun {K K'} φ => theorem_4_2_current_form_implies_naive_GC p h K K' φ)

/-- **反例が 1 組でも出れば現行形は偽**。

[8] (Jarden-Ritter, On the Characterization of Local Fields by their Absolute Galois
Groups, J. Number Th. 11 (1979), pp. 1-13) が与えるのはまさにこの仮説——絶対 Galois 群が
位相群として同型なのに ℚ_p-代数として同型でない 2 つの p 進局所体である。
その構成は本リポジトリには無いので、本補題は「その日が来たら現行形が倒れる」ことだけを
述べる(`Check/PGC/Prop12ForallRD.lean::not_forall_RD_recoverable_of_nonisomorphic` と同型)。 -/
theorem not_theorem_4_2_current_form_of_nonisomorphic {K K' : PAdicLocalField p}
    (φ : ContinuousMulEquiv K.absGal K'.absGal)
    (hne : ¬ Nonempty (K.carrier ≃ₐ[ℚ_[p]] K'.carrier)) :
    ¬ ∀ (RF : RamificationFiltration p) (hnat : IsNaturalFiltration RF)
      (K K' : PAdicLocalField p),
      Function.Bijective (naturalOuterIso RF hnat (K := K) (K' := K')) :=
  fun h => hne (theorem_4_2_current_form_implies_naive_GC p h K K' φ)

end ABC3.Check.PGC
