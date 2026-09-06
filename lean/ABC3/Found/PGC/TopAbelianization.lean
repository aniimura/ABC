import Mathlib.Topology.Algebra.Group.TopologicalAbelianization
import Mathlib.Topology.Algebra.ClopenNhdofOne
import Mathlib.Topology.Algebra.Nonarchimedean.TotallyDisconnected
import Mathlib.Topology.Algebra.ContinuousMonoidHom

/-!
# 位相的アーベル化 `Γ^{ab} = Γ ⧸ ⁅Γ,Γ⁆‾` —— 副有限性と移送

[pGC] Proposition 1.1 への**経路 Λ**(円分子 `Λ(Γ_K) ≅ μ_{p^n}` を経由する道)の
最初の節点 Λ1。`Found/PGC/TorsionCyclotome.lean`(Λ2)と
`Found/PGC/CyclotomeTransport.lean`(Λ10)が本ファイルの上に乗る。

## ★`Abelianization Γ` を使ってはいけない(設計の要)

`Abelianization Γ = Γ ⧸ commutator Γ` は**代数的**交換子群による商で、位相的閉包を取らない。
副有限群では有限指数部分群が開とは限らないので、この商を使うと

* 商が位相群として Hausdorff にならない(`commutator Γ` は閉とは限らない)、
* したがって「`Γ^{ab}` の `p^n` 捩れ」が副有限群論の言葉で扱えない、

という穴が開く。経路 C を設計したときに同じ罠を潰してある
(`ResearchPaper/pgc-goal.md` の「✗ `Abelianization Γ`」)。
本ファイルは終始 mathlib の

  `TopologicalAbelianization G = G ⧸ (commutator G).topologicalClosure`

だけを使う。

## 本ファイルが確立したもの

1. **副有限性** —— `G` が副有限(コンパクト・Hausdorff・完全不連結な位相群)なら
   `G^{ab}` も副有限。Hausdorff 性は `(commutator G).topologicalClosure` が閉であることから、
   完全不連結性は「副有限群の**閉**正規部分群による商は非アルキメデス的」
   (`instNonarchimedeanGroupQuotient`)を経由して出す。
   ★一般には完全不連結空間の商は完全不連結にならない(Cantor 集合 ↠ `[0,1]`)。
   位相群の**部分群**による商であることが効いている。
2. **移送** —— `α : Γ ≃ₜ* Γ'` は `Γ^{ab} ≃ₜ* Γ'^{ab}` を誘導する
   (`topAbelianizationCME`)。鍵は `α` が同相であることから
   `α(S‾) = α(S)‾`(`map_topologicalClosure_of_homeo`)が言えることで、
   ここでも「位相的」閉包であることが本質的に効く。

## 逸脱

無し。原典(pGC §1)は `Γ_K^{ab}` を局所類体論の文脈で暗黙に使っており、
位相的アーベル化であることは明示していないが、副有限群論の標準的な読み方である。
-/

namespace ABC3.Found.PGC

/-! ## 1. 交換子群の位相的閉包 -/

section Closure

variable (G : Type*) [Group G] [TopologicalSpace G] [IsTopologicalGroup G]

/-- `(commutator G).topologicalClosure` は閉。`T2Space (G ⧸ ·)` の合成のために
インスタンスにしておく(mathlib は `Subgroup.isClosed_topologicalClosure` を
定理としてしか持っていない)。 -/
instance instIsClosedTopologicalClosureCommutator :
    IsClosed (((commutator G).topologicalClosure : Subgroup G) : Set G) :=
  Subgroup.isClosed_topologicalClosure _

end Closure

/-! ## 2. 副有限群の閉正規部分群による商は非アルキメデス的 -/

/-- **副有限群 `G` と正規部分群 `N` について `G ⧸ N` は非アルキメデス的**、すなわち
`1` の近傍系が開部分群で細分される。

`U ∈ 𝓝 1` の引き戻しは `G` の `1` の開近傍なので、副有限性
(`ProfiniteGrp.exist_openNormalSubgroup_sub_open_nhds_of_one`)から開正規部分群
`H` が入る。商写像は開写像(`QuotientGroup.isOpenQuotientMap_mk`)なので
`H` の像は開部分群で、`U` に入る。

★`N` が閉であることはここでは要らない(Hausdorff 性のためだけに使う)。 -/
instance instNonarchimedeanGroupQuotient {G : Type*} [Group G] [TopologicalSpace G]
    [IsTopologicalGroup G] [CompactSpace G] [TotallyDisconnectedSpace G]
    (N : Subgroup G) [N.Normal] : NonarchimedeanGroup (G ⧸ N) := by
  refine ⟨fun U hU => ?_⟩
  have hpre : QuotientGroup.mk ⁻¹' U ∈ nhds (1 : G) :=
    (QuotientGroup.isOpenQuotientMap_mk (N := N)).continuous.continuousAt hU
  obtain ⟨V, hVsub, hVopen, hV1⟩ := mem_nhds_iff.mp hpre
  obtain ⟨H, hH⟩ := ProfiniteGrp.exist_openNormalSubgroup_sub_open_nhds_of_one hVopen hV1
  refine ⟨⟨(H : Subgroup G).map (QuotientGroup.mk' N), ?_⟩, ?_⟩
  · show IsOpen ((((H : Subgroup G).map (QuotientGroup.mk' N)) : Subgroup (G ⧸ N)) : Set (G ⧸ N))
    rw [Subgroup.coe_map]
    exact QuotientGroup.isOpenQuotientMap_mk.isOpenMap _ H.isOpen'
  · rintro _ ⟨x, hx, rfl⟩
    exact hVsub (hH hx)

/-! ## 3. `Γ^{ab}` の副有限性 -/

section Profinite

variable (G : Type*) [Group G] [TopologicalSpace G] [IsTopologicalGroup G]

/-- `G` がコンパクトなら `G^{ab}` もコンパクト。 -/
theorem compactSpace_topologicalAbelianization [CompactSpace G] :
    CompactSpace (TopologicalAbelianization G) := inferInstance

/-- `G` が Hausdorff なら `G^{ab}` も Hausdorff。
★これが**位相的**アーベル化を取る第一の理由:`Abelianization G` では
`commutator G` が閉とは限らないので Hausdorff にならない。 -/
theorem t2Space_topologicalAbelianization [T2Space G] :
    T2Space (TopologicalAbelianization G) := inferInstance

/-- **`G` が副有限なら `G^{ab}` も副有限**(完全不連結性)。 -/
theorem totallyDisconnectedSpace_topologicalAbelianization
    [CompactSpace G] [T2Space G] [TotallyDisconnectedSpace G] :
    TotallyDisconnectedSpace (TopologicalAbelianization G) := inferInstance

end Profinite

/-! ## 4. 移送——位相群の同型はアーベル化の同型を誘導する -/

/-- 群同型は交換子群を交換子群に写す。 -/
theorem map_commutator_of_mulEquiv {G G' : Type*} [Group G] [Group G'] (f : G ≃* G') :
    (commutator G).map f.toMonoidHom = commutator G' := by
  rw [commutator, Subgroup.map_commutator, Subgroup.map_top_of_surjective _ f.surjective]
  rfl

/-- **位相群の同型は位相的閉包と可換**——`α(S‾) = α(S)‾`。
`α` が同相であること(`Homeomorph.image_closure`)だけを使う。 -/
theorem map_topologicalClosure_of_homeo {G G' : Type*} [Group G] [Group G']
    [TopologicalSpace G] [TopologicalSpace G'] [IsTopologicalGroup G] [IsTopologicalGroup G']
    (f : G ≃ₜ* G') (S : Subgroup G) :
    (S.topologicalClosure).map f.toMulEquiv.toMonoidHom
      = (S.map f.toMulEquiv.toMonoidHom).topologicalClosure := by
  apply SetLike.coe_injective
  show ⇑f.toMulEquiv.toMonoidHom '' (closure (S : Set G))
    = closure (⇑f.toMulEquiv.toMonoidHom '' (S : Set G))
  exact f.toHomeomorph.image_closure _

/-- 位相群の同型は「交換子群の位相的閉包」を対応させる。 -/
theorem map_topCommutator {G G' : Type*} [Group G] [Group G']
    [TopologicalSpace G] [TopologicalSpace G'] [IsTopologicalGroup G] [IsTopologicalGroup G']
    (α : G ≃ₜ* G') :
    ((commutator G).topologicalClosure).map α.toMulEquiv.toMonoidHom
      = (commutator G').topologicalClosure := by
  rw [map_topologicalClosure_of_homeo, map_commutator_of_mulEquiv]

/-- **★★★★★★★★`α : Γ ≃ₜ* Γ'` が誘導する `Γ^{ab} ≃ₜ* Γ'^{ab}`**。

群同型の部分は `QuotientGroup.congr`、連続性は商写像が開商写像であること
(`QuotientGroup.isOpenQuotientMap_mk`)から両方向とも出る。 -/
noncomputable def topAbelianizationCME {G G' : Type*} [Group G] [Group G']
    [TopologicalSpace G] [TopologicalSpace G'] [IsTopologicalGroup G] [IsTopologicalGroup G']
    (α : G ≃ₜ* G') : TopologicalAbelianization G ≃ₜ* TopologicalAbelianization G' where
  toMulEquiv := QuotientGroup.congr _ _ α.toMulEquiv (map_topCommutator α)
  continuous_toFun := by
    rw [QuotientGroup.isOpenQuotientMap_mk.isQuotientMap.continuous_iff]
    exact (QuotientGroup.isOpenQuotientMap_mk.continuous).comp α.continuous_toFun
  continuous_invFun := by
    rw [QuotientGroup.isOpenQuotientMap_mk.isQuotientMap.continuous_iff]
    exact (QuotientGroup.isOpenQuotientMap_mk.continuous).comp α.continuous_invFun

@[simp] theorem topAbelianizationCME_mk {G G' : Type*} [Group G] [Group G']
    [TopologicalSpace G] [TopologicalSpace G'] [IsTopologicalGroup G] [IsTopologicalGroup G']
    (α : G ≃ₜ* G') (g : G) :
    topAbelianizationCME α (QuotientGroup.mk g) = QuotientGroup.mk (α g) := rfl

/-- 合成の関手性。 -/
theorem topAbelianizationCME_trans {G G' G'' : Type*} [Group G] [Group G'] [Group G'']
    [TopologicalSpace G] [TopologicalSpace G'] [TopologicalSpace G'']
    [IsTopologicalGroup G] [IsTopologicalGroup G'] [IsTopologicalGroup G'']
    (α : G ≃ₜ* G') (β : G' ≃ₜ* G'') (x : TopologicalAbelianization G) :
    topAbelianizationCME β (topAbelianizationCME α x) = topAbelianizationCME (α.trans β) x := by
  induction x using QuotientGroup.induction_on with
  | H g => rfl

/-- 恒等の関手性。 -/
theorem topAbelianizationCME_refl {G : Type*} [Group G] [TopologicalSpace G] [IsTopologicalGroup G]
    (x : TopologicalAbelianization G) :
    topAbelianizationCME (ContinuousMulEquiv.refl G) x = x := by
  induction x using QuotientGroup.induction_on with
  | H g => rfl

end ABC3.Found.PGC
