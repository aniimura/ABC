import ABC3.Check.PGC.Prop22Degenerate
import ABC3.Check.PGC.Prop12Degenerate
import ABC3.Found.PGC.Prop22FixedForm

/-!
# [pGC] Proposition 2.2 の**現行**形もまだ偽——自由な型族は同型不変でない

原文 (pGC p.5):

> Suppose that we are given the following group-theoretic data: the topological group Γ_K,
> together with the indexed filtration Γ_K^v for all v > 0. Then the Γ_K-modules O[scr]_K[bar], and
> K[bar]∧ can be recovered group-theoretically from this group-theoretic data.

`Check/PGC/Prop22Degenerate.lean`(2026-09-05)は、旧形の作用が公理ゼロの `SMul`
だったことを反例で示し、`DistribMulAction` へ強めて修理した。
★本ファイルはその**修理後**の形がまだ偽であることを示す。

## 反例(`prop_2_2_free_form_false`、`sorry` 無し)

今度は作用は病的でない——**自明な作用**である。落ちるのは型族の側:

* `Skeleton/PGC/Section2.lean::prop_2_2` は `O_K̄`・`K̄^` を
  **`PAdicLocalField p` の項の関数**として全称量化している。
* ところが `Γ_K ≅ Γ_K'`(位相群)なのに `K ≠ K'` である項の対が実在する
  ——`Check/PGC/Prop12Degenerate.lean` の `twistedField p` と `selfField p`
  (台の型を `ℚ_[p]` のままにして体構造だけを `x ↦ -x` に沿って捻ったもの)。
* そこで「`K = twistedField p` のときだけ `ℤ`、他は `0`」という型族を取る
  (`Found/PGC/Prop22FixedForm.lean::condSubgroup`)。作用は自明な作用でよいので
  `DistribMulAction` の公理はすべて満たされる。
* `α := twistedGalEquiv p` を当てると `ℤ ≃+ 0` が要求されて落ちる。

★落ちているのは原典の数学ではなく、**我々の形式化**である:

> 型族 `IntKbar` は `PAdicLocalField p` の項の関数であって、同型類の関数ではない。
> 同型不変性を課していないので、「Γ_K から回復できる」は最初から成り立ちようがない。

すなわち **D13**(`Skeleton/PGC/Section1.lean` の Prop 1.2 で `∀ RD` をやめて
実物に固定した判断)と**同じ退化**であり、落とした条件も同じ「同型不変性」である。

## ★修理は「実物に固定する」

D13 と同じく、修理は**自由な型族をやめて実物に固定する**ことである。実物は
`Found/PGC/AbsClosureModules.lean` に構成済みで(`IntKbar K = ↥(absClosureInt K)`、
`CompKbar K = closureCompletion K`)、固定した形の還元は
`Found/PGC/Prop22FixedForm.lean::prop_2_2_real_of_isometric` に在る。

★★**修理が「弱めただけ」ではないことの証拠**:
`intKbar_transport_twisted` / `compKbar_transport_twisted`(本ファイル末尾)は、
★**上の反例が使うのとまったく同じ α**(`twistedGalEquiv p`)について、
実物 `𝒪_{K̄}`・`ℂ_K` の同変加法同型が**無条件に存在する**ことを示す。
すなわち固定版は反例の場所を避けているのではなく、そこで実際に成り立っている。

★これで「落とした条件は、主張を偽にするか自明にするかのどちらかになる」例は
`InertiaDegeneracy`・`Theorem42Degenerate`・`Def32Degenerate`・`Cor33Degenerate`・
`Prop22Degenerate`・`Prop12Degenerate`・`Prop12ForallRD`・本件 で 8 つ目。

**これは原典の主張ではない**(我々のモデルと器具についての事実)ので `.src` を持たない。
-/

namespace ABC3.Check.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC ABC3.Found.PGC

/-- **★★★★★★[pGC] Proposition 2.2 の現行形も偽**——作用は `DistribMulAction`
(自明な作用)で公理をすべて満たすが、自由な型族が同型不変でないので落ちる。

★`Check/PGC/Prop22Degenerate.lean::prop_2_2_statement_false` とは**別の**反例である。
あちらは作用の病理(非可換性を使う)、こちらは型族の病理(項の非等号を使う)。 -/
theorem prop_2_2_free_form_false (p : ℕ) [Fact p.Prime] :
    ¬ (∀ (_RF : RamificationFiltration p)
        (IntKbar CompKbar : PAdicLocalField p → Type)
        (i1 : ∀ K, AddCommGroup (IntKbar K))
        (i2 : ∀ K, DistribMulAction K.absGal (IntKbar K))
        (i3 : ∀ K, AddCommGroup (CompKbar K))
        (i4 : ∀ K, DistribMulAction K.absGal (CompKbar K)),
        @RecoverableAsAddModule p _ IntKbar i1 i2
          ∧ @RecoverableAsAddModule p _ CompKbar i3 i4) := by
  intro h
  obtain ⟨key, -⟩ := h (topFiltration' p)
      (fun K => ↥(condSubgroup (K = twistedField p)))
      (fun K => ↥(condSubgroup (K = twistedField p)))
      (fun _ => inferInstance) (fun _ => constDistribMulAction _ _)
      (fun _ => inferInstance) (fun _ => constDistribMulAction _ _)
  obtain ⟨φ, -⟩ := key (K := twistedField p) (K' := selfField p) (twistedGalEquiv p)
  exact false_of_addEquiv_condSubgroup rfl (fun hc => twistedField_ne_selfField p hc.symm) φ

/-! ## ★修理後は、反例と同じ α の上で実際に成り立つ -/

/-- ★★**自由版が落ちるその α において、実物 `𝒪_{K̄}` は移送できる**(無条件)。

★α は `twistedGalEquiv p` そのものである——`Check/PGC/Prop12Degenerate.lean` の
定義がまさに `galContinuousMulEquiv (twistedAlgEquiv p)` だからである。
ここで `twistedGalEquiv p` と書かずに右辺を展開して書いているのは配管上の理由で、
`TwistedQp p` が `ℚ_[p]` の型同義語であるために単一化が `whnf` で止まるのを避けるため
(`tools/lean-idioms.md` #262)。★`(K := …) (K' := …)` を明示すれば止まらない。 -/
theorem intKbar_transport_twisted (p : ℕ) [Fact p.Prime] :
    ∃ φ : IntKbar (twistedField p) ≃+ IntKbar (selfField p),
      ∀ (g : (twistedField p).absGal) (x : IntKbar (twistedField p)),
        φ (g • x)
          = ((galContinuousMulEquiv (K := twistedField p) (K' := selfField p)
              (twistedAlgEquiv p)).toMulEquiv g) • φ x :=
  intKbar_transport_galContinuousMulEquiv (K := twistedField p) (K' := selfField p)
    (twistedAlgEquiv p)

/-- ★★同じ α について `ℂ_K` も移送できる(無条件)。 -/
theorem compKbar_transport_twisted (p : ℕ) [Fact p.Prime] :
    ∃ φ : CompKbar (twistedField p) ≃+ CompKbar (selfField p),
      ∀ (g : (twistedField p).absGal) (x : CompKbar (twistedField p)),
        φ (g • x)
          = ((galContinuousMulEquiv (K := twistedField p) (K' := selfField p)
              (twistedAlgEquiv p)).toMulEquiv g) • φ x :=
  compKbar_transport_galContinuousMulEquiv (K := twistedField p) (K' := selfField p)
    (twistedAlgEquiv p)

end ABC3.Check.PGC
