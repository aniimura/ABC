import ABC3.Found.PGC.ReciprocityAlphaTransport
import ABC3.Found.PGC.LocalClassFieldTheory
import ABC3.Found.PGC.TopAbelianization

/-!
# N1 —— `ker(Art_K) ⊓ I_K` は群論的に標準である

## ★★★0. 結論を先に —— 4 行

1. ★**N1 の字面は真**である。しかも右辺の `⊓ I_K` は**要らない**:

   ```
   ker(Art_π) ⊓ I_K  =  Gal(K̄/K^ab)  =  ‾⁅Γ_K, Γ_K⁆
   ```

   (`ker_reciprocityUnits_inf_absInertia` /
   `ker_reciprocityUnits_inf_absInertia_eq_topCommutator`)。
   ★右辺には **`π` も `f` も局所体の語彙も出てこない** —— `Γ_K` の**群と位相だけ**である。
2. ★★★**したがって `ArtinKerTransport` の惰性版は仮説なしで成り立つ**
   (`artinKerInertiaTransport`)。★α には「位相群同型である」以外の仮定が要らない。
3. ★★★★★**`ReciprocityAlphaTransport.lean` の結論の鎖は、惰性群に制限した形で
   丸ごと無条件になる** —— `unitsTransportInertia` / `map_principalUnits_unitsTransportInertia` /
   `integersTransportInertia : 𝒪[K] ≃+ 𝒪[K']`。
   ★仮説は「`α` が濾過つき位相群同型であること」だけになった。
4. ★逆に、`ker(Art_K)` を `I_K` に制限しない版が偽であることは変わらない
   (`ReciprocityAlphaTransport.lean::artinKerTransport_refl_iff_lubinTateClosure_eq`)。
   ★★制限が本質的である。

## ★★1. なぜ真か —— 3 行の計算

Galois 対応で全部書き換えるだけである:

```
ker(Art_π) ⊓ I_K
  = Gal(K̄/K_π) ⊓ Gal(K̄/K^ur)        (ker_reciprocityUnits / absInertia の定義)
  = Gal(K̄/(K_π ⊔ K^ur))              (IntermediateField.fixingSubgroup_sup、mathlib)
  = Gal(K̄/K^ab)                      (★局所 Kronecker-Weber)
  = ‾⁅Γ_K, Γ_K⁆                       (K^ab の定義 + 閉部分群の Galois 対応)
```

★★`π` 依存性は 3 行目で消える。`K_π` は `π` に依るが `K_π ⊔ K^ur` は依らない。
★これが `ReciprocityAlphaTransport` の反例
(`K = ℚ_p`、`π = p` vs `π′ = -p`、`K_p ≠ K_{-p}`)を回避する仕組みである ——
`K_p ⊔ K^ur = K_{-p} ⊔ K^ur = ℚ_p^ab` だからである。

★`⊓ I_K` を落とせるのは `K^ur ≤ K^ab` (`unramifiedClosure_le_abelianClosure`)、
すなわち `Gal(K̄/K^ab) ≤ I_K` だからである。★持ち場の字面
「`ker(Art_K) ∩ I_K = Gal(K̄/K^ab) ∩ I_K`」は真だが右辺の `∩ I_K` が冗長である。

## ★★2. 在庫調査(コマンドを残す)

```
grep -nE "topologicalClosure" .cache/decl-index.txt | grep -E "commutator"
grep -nE "IntermediateField.fixingSubgroup_sup" .cache/mathlib-index.txt
grep -nE "absInertia" .cache/decl-index.txt
grep -riln "kronecker" lean/ABC3/ --include=*.lean
```

★★「無いと思ったが在った」が 4 件(すべて書く前に測って回避した):

| 書こうとしたもの | 実際に在った場所 |
|---|---|
| `map α (‾⁅G,G⁆) = ‾⁅G',G'⁆` | ★`ABC3.Found.PGC.map_topCommutator`(`TopAbelianization.lean:135`) |
| `map f ⁅G,G⁆ = ⁅G',G'⁆`(全射) | `ABC3.Found.PGC.map_commutator_of_mulEquiv`(同 `:117`) |
| `map α (closure H) = closure (map α H)` | `ABC3.Found.PGC.map_topologicalClosure_of_homeo`(同 `:124`) |
| `(K ⊔ L).fixingSubgroup = K.fS ⊓ L.fS` | `IntermediateField.fixingSubgroup_sup`(mathlib、★`FiniteDimensional` 不要) |

★局所 Kronecker-Weber も木に在った:
`abelianClosure_eq_lubinTateClosure_sup_unramifiedClosure`
(`Found/PGC/LocalClassFieldTheory.lean:723`、Yoshida08 Theorem 6.15、`sorry` 無し)。
★これが「`K_π ⊔ K^ur = K^ab`」そのものである(持ち場が別々に数えていた 2 件は同じ 1 件)。

★「在るが形が違う」も 2 件: `IntermediateField.le_iff_le` と
`instIsClosedTopologicalClosureCommutator` は索引の行に出ない明示引数を持つ
(`lean-idioms.md` #297)。

## ★★3. 抽象核(§1)—— 型クラスは `Group` 4 つだけ

`Intertwines` の部分群版である。`ReciprocityAlphaTransport.lean::transportOfKerEq`
(核が対応すれば運べる)を、部分群 `N` に制限した準同型に代入する:

* `restrictEquiv` —— `α` を `N` に制限した同型 `N ≃* N'`
* `surjective_comp_subtype` —— 像が `⊤` なら制限は全射
* `map_restrictEquiv_ker` —— `ker ⊓ N` が対応 ⟹ 制限した核が対応
* `transportOfKerInfEq` / `transportOfKerInfEq_apply` —— ★到達点
* `map_map_transportOfKerInfEq` —— ★`S ≤ N` なら部分群の族はそのまま運ばれる

★分岐・付値・Galois・Lubin-Tate・局所体は §1 に 1 語も出てこない。

## ★★4. 退化の自己検査

0. ★★**両辺が潰れていない**: `topCommutator_ne_absInertia` が
   `‾⁅Γ_K,Γ_K⁆ ≠ I_K` を示す。★これが無いと N1 は「`I_K = I_K`」という空虚な等式で
   ありうる。★中身は `I_K / ‾⁅Γ_K,Γ_K⁆ ≅ 𝒪_K^×` である。
1. ★**非空虚**: `nonempty_artinInertiaDatum` が仮定ゼロで
   `ArtinInertiaDatum K` を作る(N1 をそこに詰めている)。
2. ★★**内容がある**: `artinKerInertiaTransport` は
   `ReciprocityAlphaTransport::ArtinKerTransport` と違い恒等 `α` に限らず
   任意の位相群同型で成り立つ。★しかも `π`・`π′` が違ってよい
   (`ker_inf_absInertia_eq_of_lubinTate` —— 制限しない版が偽になるまさにその状況で真)。
3. ★★**制限は落とせない**: `⊓ I_K` を外すと偽である
   (`ReciprocityAlphaTransport::artinKerTransport_refl_iff_lubinTateClosure_eq` +
   `LubinTateUniformizerIndependence.lean` の反例)。★本ファイルはその反例を
   回避したのではなく、主張を反例が通らない形に絞ったのである。
4. ★★**まだ Proposition 2.2 全体ではない**: 得られるのは底体の `𝒪_K ≃+ 𝒪_{K'}` で、
   `𝒪_{K̄}` の Γ-同変同型ではない(`ReciprocityAlphaTransport` §4-4 と同じ限界)。

## 逸脱の記録

1. 持ち場の字面は `ker(Art_K) ∩ I_K = Gal(K̄/K^ab) ∩ I_K` だったが、
   右辺の `∩ I_K` は冗長なので落とした(`K^ur ≤ K^ab` より
   `Gal(K̄/K^ab) ≤ I_K`)。★冗長版も `ker_reciprocityUnits_inf_absInertia_inf` として残す。
2. 原典 pGC は `Art` を `Γ^ab_K` の言葉で述べるが、木の `reciprocityUnits` は
   `𝒪_K^×` 成分だけを取る(`ReciprocityAlphaTransport` の逸脱 2 を踏襲)。
3. `Skeleton/**` は 1 行も触っていない。`Found/PGC/LubinTate*.lean` も読むだけ。
   `ReciprocityAlphaTransport.lean` も読むだけ(import して使うのみ)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open IsLocalRing
open scoped NNReal Valued

/-! ## §1 ★★★抽象核 —— 純群論。分岐・付値・Galois・局所体は 1 語も出ない -/

section RestrictCore

variable {G G' H H' : Type*} [Group G] [Group G'] [Group H] [Group H']

/-- `α : G ≃* G'` を部分群 `N` に制限した同型 `N ≃* N'`。 -/
def restrictEquiv (α : G ≃* G') {N : Subgroup G} {N' : Subgroup G'}
    (hα : Subgroup.map (α : G →* G') N = N') : N ≃* N' :=
  (α.subgroupMap N).trans (MulEquiv.subgroupCongr hα)

@[simp] theorem coe_restrictEquiv (α : G ≃* G') {N : Subgroup G} {N' : Subgroup G'}
    (hα : Subgroup.map (α : G →* G') N = N') (x : N) :
    ((restrictEquiv α hα x : N') : G') = α (x : G) := rfl

/-- ★像が `⊤` なら、部分群への制限は全射。 -/
theorem surjective_comp_subtype {φ : G →* H} {N : Subgroup G} (hN : Subgroup.map φ N = ⊤) :
    Function.Surjective (φ.comp N.subtype) := by
  intro h
  have hmem : h ∈ Subgroup.map φ N := hN ▸ Subgroup.mem_top h
  obtain ⟨g, hg, rfl⟩ := hmem
  exact ⟨⟨g, hg⟩, rfl⟩

/-- ★★`ker φ ⊓ N` が対応するなら、`N` に制限した準同型の核が対応する。 -/
theorem map_restrictEquiv_ker (φ : G →* H) (φ' : G' →* H') (α : G ≃* G')
    {N : Subgroup G} {N' : Subgroup G'} (hα : Subgroup.map (α : G →* G') N = N')
    (hker : Subgroup.map (α : G →* G') (φ.ker ⊓ N) = φ'.ker ⊓ N') :
    Subgroup.map ((restrictEquiv α hα : N ≃* N') : N →* N') (φ.comp N.subtype).ker
      = (φ'.comp N'.subtype).ker := by
  ext y
  simp only [Subgroup.mem_map, MonoidHom.mem_ker, MonoidHom.coe_comp, Function.comp_apply,
    Subgroup.coe_subtype, MonoidHom.coe_coe]
  constructor
  · rintro ⟨x, hx, rfl⟩
    have hm : (x : G) ∈ φ.ker ⊓ N := ⟨hx, x.2⟩
    have hmm : α (x : G) ∈ φ'.ker ⊓ N' := hker ▸ ⟨_, hm, rfl⟩
    simpa [coe_restrictEquiv] using hmm.1
  · intro hy
    have hm : (y : G') ∈ φ'.ker ⊓ N' := ⟨hy, y.2⟩
    rw [← hker] at hm
    obtain ⟨g, hg, hgy⟩ := hm
    refine ⟨⟨g, hg.2⟩, hg.1, ?_⟩
    exact Subtype.ext (by simpa [coe_restrictEquiv] using hgy)

variable (φ : G →* H) (φ' : G' →* H') {N : Subgroup G} {N' : Subgroup G'}
  (hN : Subgroup.map φ N = ⊤) (hN' : Subgroup.map φ' N' = ⊤)
  (α : G ≃* G') (hα : Subgroup.map (α : G →* G') N = N')
  (hker : Subgroup.map (α : G →* G') (φ.ker ⊓ N) = φ'.ker ⊓ N')

/-- ★★★**抽象核の到達点** —— `N` の上で全射で、`ker ⊓ N` が対応すれば運べる。

★`ReciprocityAlphaTransport.lean::transportOfKerEq`(核が対応すれば運べる)を
`φ ∘ N.subtype` に代入するだけである。 -/
noncomputable def transportOfKerInfEq : H ≃* H' :=
  transportOfKerEq (φ.comp N.subtype) (surjective_comp_subtype hN)
    (φ'.comp N'.subtype) (surjective_comp_subtype hN')
    (restrictEquiv α hα) (map_restrictEquiv_ker φ φ' α hα hker)

/-- ★★α-同変性そのもの(ただし `N` の上でだけ)。 -/
theorem transportOfKerInfEq_apply {g : G} (hg : g ∈ N) :
    transportOfKerInfEq φ φ' hN hN' α hα hker (φ g) = φ' (α g) :=
  transportOfKerEq_intertwines (φ.comp N.subtype) (surjective_comp_subtype hN)
    (φ'.comp N'.subtype) (surjective_comp_subtype hN')
    (restrictEquiv α hα) (map_restrictEquiv_ker φ φ' α hα hker) ⟨g, hg⟩

/-- ★★`N` の内側の部分群の族はそのまま運ばれる。 -/
theorem map_map_transportOfKerInfEq {S : Subgroup G} (hS : S ≤ N) :
    Subgroup.map ((transportOfKerInfEq φ φ' hN hN' α hα hker : H ≃* H') : H →* H')
        (Subgroup.map φ S)
      = Subgroup.map φ' (Subgroup.map (α : G →* G') S) := by
  ext z
  simp only [Subgroup.mem_map, MonoidHom.coe_coe, exists_exists_and_eq_and]
  constructor
  · rintro ⟨g, hg, rfl⟩
    exact ⟨g, hg, (transportOfKerInfEq_apply φ φ' hN hN' α hα hker (hS hg)).symm⟩
  · rintro ⟨g, hg, rfl⟩
    exact ⟨g, hg, transportOfKerInfEq_apply φ φ' hN hN' α hα hker (hS hg)⟩

end RestrictCore

/-! ## §1b 準抽象核 —— 閉部分群の Galois 対応(局所体は出てこない) -/

section GaloisCore

variable {k E : Type*} [Field k] [Field E] [Algebra k E]

/-- ★閉部分群なら `fixingSubgroup ∘ fixedField` は恒等。

`≤` は `AbelianClosureSplit.lean::fixingSubgroup_fixedField_le_topologicalClosure` +
閉包の最小性、`≥` は Galois 接続 `IntermediateField.le_iff_le`。
★mathlib の `InfiniteGalois.fixingSubgroup_fixedField` は `ClosedSubgroup` の
束ね直しを要求するので、`IsClosed` を直接受け取る形にした。 -/
theorem fixingSubgroup_fixedField_of_isClosed [IsGalois k E]
    {H : Subgroup (E ≃ₐ[k] E)} (hH : IsClosed ((H : Subgroup (E ≃ₐ[k] E)) : Set (E ≃ₐ[k] E))) :
    (IntermediateField.fixedField H).fixingSubgroup = H :=
  le_antisymm
    ((fixingSubgroup_fixedField_le_topologicalClosure H).trans
      (Subgroup.topologicalClosure_minimal H le_rfl hH))
    ((IntermediateField.le_iff_le H (IntermediateField.fixedField H)).mp le_rfl)

end GaloisCore

variable {p : ℕ} [Fact p.Prime]

/-- ★★`Gal(K̄/K^ab) = ‾⁅Γ_K, Γ_K⁆` —— `K^ab` の側から `Γ_K` の側へ戻る。 -/
theorem fixingSubgroup_abelianClosure (K : PAdicLocalField p) :
    (abelianClosure K).fixingSubgroup = (commutator K.absGal).topologicalClosure := by
  haveI := isGalois_closure K
  rw [abelianClosure_def]
  exact fixingSubgroup_fixedField_of_isClosed (Subgroup.isClosed_topologicalClosure _)

def fixingSubgroup_abelianClosure.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Theorem 6.15", sectionId := "thm-6-15" }

/-- ★`Gal(K̄/K^ab) ≤ I_K`(`K^ur ≤ K^ab` の Galois 対応)。 -/
theorem fixingSubgroup_abelianClosure_le_absInertia (K : PAdicLocalField p) :
    (abelianClosure K).fixingSubgroup ≤ absInertia K :=
  IntermediateField.fixingSubgroup_le (unramifiedClosure_le_abelianClosure K)

/-! ## §2 ★★★★★N1 本体 -/

section N1

variable (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    [Normal K.carrier (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf)]

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**N1** —— `ker(Art_π) ⊓ I_K = Gal(K̄/K^ab)`。

★★右辺に `π` も `f` も出てこない。左辺は Lubin-Tate の選択に依るように見えるが、
`⊓ I_K` を取ると依存性が消える。

証明は 4 つの書き換えだけ:
`ker Art_π = Gal(K̄/K_π)`(`ker_reciprocityUnits`)、
`I_K = Gal(K̄/K^ur)`(`absInertia` の定義)、
`Gal(K̄/A) ⊓ Gal(K̄/B) = Gal(K̄/(A ⊔ B))`(`IntermediateField.fixingSubgroup_sup`)、
`K_π ⊔ K^ur = K^ab`(★局所 Kronecker-Weber)。 -/
theorem ker_reciprocityUnits_inf_absInertia :
    (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf).ker ⊓ absInertia K
      = (abelianClosure K).fixingSubgroup := by
  rw [ker_reciprocityUnits, show absInertia K = (unramifiedClosure K).fixingSubgroup from rfl,
    ← IntermediateField.fixingSubgroup_sup,
    abelianClosure_eq_lubinTateClosure_sup_unramifiedClosure K hq hπmax hπne0 f hf0 hf1 hf]

def ker_reciprocityUnits_inf_absInertia.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- ★★★★N1 の群論的な形 —— 右辺は `Γ_K` の群と位相だけで書けている。

```
ker(Art_π) ⊓ I_K = ‾⁅Γ_K, Γ_K⁆
```

★これが「群論的に標準」の中身である。★体も付値も分岐も右辺に現れない。 -/
theorem ker_reciprocityUnits_inf_absInertia_eq_topCommutator :
    (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf).ker ⊓ absInertia K
      = (commutator K.absGal).topologicalClosure := by
  rw [ker_reciprocityUnits_inf_absInertia K hq hπmax hπne0 f hf0 hf1 hf,
    fixingSubgroup_abelianClosure K]

def ker_reciprocityUnits_inf_absInertia_eq_topCommutator.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- ★持ち場の字面そのまま(右辺の `⊓ I_K` つき)。★冗長だが、消費側が探しに来る形。 -/
theorem ker_reciprocityUnits_inf_absInertia_inf :
    (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf).ker ⊓ absInertia K
      = (abelianClosure K).fixingSubgroup ⊓ absInertia K := by
  rw [ker_reciprocityUnits_inf_absInertia K hq hπmax hπne0 f hf0 hf1 hf,
    inf_eq_left.mpr (fixingSubgroup_abelianClosure_le_absInertia K)]

end N1

/-! ## §3 ★★★★★仮説なしの α-移送 -/

/-- ★★**惰性版 Artin データ** —— `ArtinFilteredDatum` に
「`ker ⊓ I_K` が群論的に標準である」を足しただけ。

★この 1 行が `ReciprocityAlphaTransport.lean` の `ArtinKerTransport` を**不要にする**。 -/
structure ArtinInertiaDatum (K : PAdicLocalField p) extends ArtinFilteredDatum K where
  /-- ★N1: 核と惰性群の交わりは閉交換子群である。 -/
  ker_inf_absInertia :
    toArtinFilteredDatum.art.ker ⊓ absInertia K = (commutator K.absGal).topologicalClosure

def ArtinInertiaDatum.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- ★★★★**惰性版 Artin データは仮定ゼロで存在する**。

`nonempty_artinFilteredDatum` と同じ組み上げに、N1 を 1 つ足すだけ。 -/
theorem nonempty_artinInertiaDatum (K : PAdicLocalField p) : Nonempty (ArtinInertiaDatum K) := by
  haveI := isAdicComplete_valuationRing K
  haveI := valuationRing_isDVR K
  obtain ⟨ϖ, hϖirr⟩ := IsDiscreteValuationRing.exists_irreducible (𝒪[K.carrier])
  have hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {ϖ} :=
    (IsDiscreteValuationRing.irreducible_iff_uniformizer ϖ).mp hϖirr
  have hπne0 : ϖ ≠ 0 := hϖirr.ne_zero
  have hq : Fintype.card 𝓀[K.carrier] = p ^ (absoluteInertiaDegree K) := by
    rw [← Nat.card_eq_fintype_card]
    exact residueCard_eq_pow K
  obtain ⟨f, hf0, hf1, hf⟩ := exists_lubinTateSeries (A := 𝒪[K.carrier]) hq hπmax
  haveI : Normal K.carrier (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf) :=
    normal_lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf
  exact ⟨{ toArtinFilteredDatum := lubinTateArtinFilteredDatum K hq hπmax hπne0 f hf0 hf1 hf
           ker_inf_absInertia :=
             ker_reciprocityUnits_inf_absInertia_eq_topCommutator K hq hπmax hπne0 f hf0 hf1 hf }⟩

def nonempty_artinInertiaDatum.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**核の対応は仮説なしで成り立つ**(惰性群に制限した版)。

```
α (ker Art_K ⊓ I_K) = ker Art_{K'} ⊓ I_{K'}
```

★★`α` に要求するのは**位相群同型であること**だけ。★濾過つきである必要すらない。
★★これが `ReciprocityAlphaTransport.lean::ArtinKerTransport`(仮説)との違いである ——
制限しない版は `α = id` でも偽になりうるが、制限した版は**常に真**。

証明は N1 で両辺を `‾⁅Γ,Γ⁆` に書き換えて `map_topCommutator`
(`TopAbelianization.lean:135`、既に木に在った)を当てるだけ。 -/
theorem artinKerInertiaTransport {K K' : PAdicLocalField p} (D : ArtinInertiaDatum K)
    (D' : ArtinInertiaDatum K') (α : ContinuousMulEquiv K.absGal K'.absGal) :
    Subgroup.map (α.toMulEquiv : K.absGal →* K'.absGal) (D.art.ker ⊓ absInertia K)
      = D'.art.ker ⊓ absInertia K' := by
  rw [D.ker_inf_absInertia, D'.ker_inf_absInertia]
  exact map_topCommutator α

def artinKerInertiaTransport.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- ★★**濾過つき同型は惰性群を惰性群へ写す**(`v = 0` を入れるだけ)。 -/
theorem map_absInertia_filteredIso {K K' : PAdicLocalField p}
    (α : FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K')) :
    Subgroup.map ((filteredIsoEquiv α).toMulEquiv : K.absGal →* K'.absGal) (absInertia K)
      = absInertia K' := by
  have h := α.map_Gv (0 : ℝ)
  rwa [show (pgcFilteredGroup K).Gv (0 : ℝ) = (ramificationFiltration p).Gv K 0 from rfl,
    show (pgcFilteredGroup K').Gv (0 : ℝ) = (ramificationFiltration p).Gv K' 0 from rfl,
    ramificationFiltration_Gv_zero, ramificationFiltration_Gv_zero] at h

/-! ## §4 ★★★★★帰結 —— `𝒪_K ≃+ 𝒪_{K'}` が**仮説なしで**出る -/

section Consequences

variable {K K' : PAdicLocalField p} (D : ArtinInertiaDatum K) (D' : ArtinInertiaDatum K')
    (α : FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K'))

/-- ★★★**単数群の移送** `𝒪_K^× ≃* 𝒪_{K'}^×` —— ★仮説なし。

`ReciprocityAlphaTransport::unitsTransport` は `ArtinKerTransport` を要求したが、
こちらは要求しない。 -/
noncomputable def unitsTransportInertia : (𝒪[K.carrier])ˣ ≃* (𝒪[K'.carrier])ˣ :=
  transportOfKerInfEq D.art D'.art D.toArtinFilteredDatum.map_absInertia
    D'.toArtinFilteredDatum.map_absInertia (filteredIsoEquiv α).toMulEquiv
    (map_absInertia_filteredIso α) (artinKerInertiaTransport D D' (filteredIsoEquiv α))

def unitsTransportInertia.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- ★★**α-同変性**(惰性群の上で)。 -/
theorem unitsTransportInertia_art {g : K.absGal} (hg : g ∈ absInertia K) :
    unitsTransportInertia D D' α (D.art g) = D'.art (filteredIsoEquiv α g) :=
  transportOfKerInfEq_apply D.art D'.art D.toArtinFilteredDatum.map_absInertia
    D'.toArtinFilteredDatum.map_absInertia (filteredIsoEquiv α).toMulEquiv
    (map_absInertia_filteredIso α) (artinKerInertiaTransport D D' (filteredIsoEquiv α)) hg

/-- ★★★**α-同変な単数群同型は仮説なしで存在する**。

★`ReciprocityAlphaTransport::exists_unitsEquiv_iff_artinKerTransport` は
「存在 ⟺ 核が対応」という**同値**しか言えなかった(そして核の対応は偽になりうる)。
★★惰性群に制限すると、右辺が常に成り立つので**存在が無条件になる**。 -/
theorem exists_unitsEquiv_inertia :
    ∃ e : (𝒪[K.carrier])ˣ ≃* (𝒪[K'.carrier])ˣ,
      ∀ g ∈ absInertia K, e (D.art g) = D'.art (filteredIsoEquiv α g) :=
  ⟨unitsTransportInertia D D' α, fun _ hg => unitsTransportInertia_art D D' α hg⟩

def exists_unitsEquiv_inertia.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- ★運び方は一意(`Art(I_K) = 𝒪_K^×` だから)。 -/
theorem unitsTransportInertia_unique (e : (𝒪[K.carrier])ˣ ≃* (𝒪[K'.carrier])ˣ)
    (he : ∀ g ∈ absInertia K, e (D.art g) = D'.art (filteredIsoEquiv α g)) :
    e = unitsTransportInertia D D' α := by
  ext u
  obtain ⟨g, hg, rfl⟩ : ∃ g ∈ absInertia K, D.art g = u := by
    have hu : u ∈ Subgroup.map D.art (absInertia K) :=
      D.toArtinFilteredDatum.map_absInertia ▸ Subgroup.mem_top u
    simpa using hu
  rw [he g hg, unitsTransportInertia_art D D' α hg]

/-- ★★★**`U^n_K ↦ U^n_{K'}`(すべての `n`)** —— ★仮説なし。

`Γ^n_K ≤ I_K`(`n ≥ 0`)なので、惰性群に制限しても分岐濾過は全部見えている。 -/
theorem map_principalUnits_unitsTransportInertia (n : ℕ) :
    Subgroup.map ((unitsTransportInertia D D' α : (𝒪[K.carrier])ˣ ≃* (𝒪[K'.carrier])ˣ) :
        (𝒪[K.carrier])ˣ →* (𝒪[K'.carrier])ˣ)
        (principalUnits K D.unif n) = principalUnits K' D'.unif n := by
  have hle : (ramificationFiltration p).Gv K ((n : ℕ) : ℝ) ≤ absInertia K :=
    ramificationFiltration_Gv_le_absInertia K (by positivity)
  have hmap := map_map_transportOfKerInfEq D.art D'.art D.toArtinFilteredDatum.map_absInertia
    D'.toArtinFilteredDatum.map_absInertia (filteredIsoEquiv α).toMulEquiv
    (map_absInertia_filteredIso α) (artinKerInertiaTransport D D' (filteredIsoEquiv α)) hle
  rw [D.toArtinFilteredDatum.map_Gv n] at hmap
  rw [unitsTransportInertia, hmap, show Subgroup.map
      ((filteredIsoEquiv α).toMulEquiv : K.absGal →* K'.absGal)
      ((ramificationFiltration p).Gv K ((n : ℕ) : ℝ))
      = (ramificationFiltration p).Gv K' ((n : ℕ) : ℝ) from α.map_Gv ((n : ℕ) : ℝ)]
  exact D'.toArtinFilteredDatum.map_Gv n

def map_principalUnits_unitsTransportInertia.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- ★★系: `(𝒪_K/π^n)^× ≅ (𝒪_{K'}/π'^n)^×` —— ★仮説なし。 -/
noncomputable def unitsQuotientTransportInertia (n : ℕ) (hn : 1 ≤ n) :
    (𝒪[K.carrier] ⧸ Ideal.span ({D.unif ^ n} : Set 𝒪[K.carrier]))ˣ
      ≃* (𝒪[K'.carrier] ⧸ Ideal.span ({D'.unif ^ n} : Set 𝒪[K'.carrier]))ˣ :=
  ((principalUnitsQuotientEquiv K D.unif_max n hn).symm.trans
      (QuotientGroup.congr (principalUnits K D.unif n) (principalUnits K' D'.unif n)
        (unitsTransportInertia D D' α) (map_principalUnits_unitsTransportInertia D D' α n))).trans
    (principalUnitsQuotientEquiv K' D'.unif_max n hn)

/-- ★★退化していないことの witness(`n = 1`)—— 剰余体の単数群が対応する。 -/
noncomputable def residueUnitsTransportInertia :
    (𝒪[K.carrier] ⧸ Ideal.span ({D.unif ^ 1} : Set 𝒪[K.carrier]))ˣ
      ≃* (𝒪[K'.carrier] ⧸ Ideal.span ({D'.unif ^ 1} : Set 𝒪[K'.carrier]))ˣ :=
  unitsQuotientTransportInertia D D' α 1 le_rfl

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**pGC Proposition 2.2 の第一段が、仮説なしで α に沿って運べる。**

```
𝒪_K  ≃(log)  U^{2 e_K e_{K'}}_K  ≃(Art^{-1}|_{I_K})  Γ_K^{2 e_K e_{K'}}
     ≃(α)    Γ_{K'}^{2 e_K e_{K'}}  ≃(Art|_{I_{K'}})  U^{2 e_K e_{K'}}_{K'}  ≃(log)  𝒪_{K'}
```

★★`ReciprocityAlphaTransport::integersTransport` は仮説 `ArtinKerTransport` を
受け取っていた。★★**本定理はそれを受け取らない。** 使うのは
「`α` が濾過つき位相群同型である」ことだけである。 -/
noncomputable def integersTransportInertia : 𝒪[K.carrier] ≃+ 𝒪[K'.carrier] :=
  AddEquiv.toMultiplicative.symm
    ((padicLogPrincipalUnitsEquiv K D.unif_max D.unif_ne_zero
        (two_le_two_mul_absoluteRamificationIndex K')).symm.trans
      (((unitsTransportInertia D D' α).subgroupMap
          (principalUnits K D.unif
            (2 * absoluteRamificationIndex K' * absoluteRamificationIndex K))).trans
        ((MulEquiv.subgroupCongr (map_principalUnits_unitsTransportInertia D D' α _)).trans
          ((MulEquiv.subgroupCongr
              (congrArg (principalUnits K' D'.unif) (by ring))).trans
            (padicLogPrincipalUnitsEquiv K' D'.unif_max D'.unif_ne_zero
              (two_le_two_mul_absoluteRamificationIndex K))))))

def integersTransportInertia.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

end Consequences

/-- ★★★★★**仮説を 1 つも受け取らない形** —— pGC Proposition 2.2 の第一段。

濾過つき位相群同型 `α : Γ_K ≃ Γ_{K'}` があれば、`𝒪_K ≃+ 𝒪_{K'}` が在る。
★Artin データも Lubin-Tate の選択も `∃` の内側に閉じ込めてある。 -/
theorem nonempty_integers_addEquiv_of_filteredIso {K K' : PAdicLocalField p}
    (α : FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K')) :
    Nonempty (𝒪[K.carrier] ≃+ 𝒪[K'.carrier]) := by
  obtain ⟨D⟩ := nonempty_artinInertiaDatum K
  obtain ⟨D'⟩ := nonempty_artinInertiaDatum K'
  exact ⟨integersTransportInertia D D' α⟩

def nonempty_integers_addEquiv_of_filteredIso.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-! ## §5 ★★退化していないことの検査 + 制限しない版との対比 -/

/-- `𝒪_K^×` は自明でない(`-1 ≠ 1`。`𝒪_K` は標数 0)。
★mathlib に `Nontrivial (𝒪[K])ˣ` のインスタンスは**無い**
(`infer_instance` が `failed to synthesize` を返す)ので、ここで作る。 -/
theorem nontrivial_units (K : PAdicLocalField p) : Nontrivial (𝒪[K.carrier])ˣ := by
  refine ⟨⟨1, -1, fun h => ?_⟩⟩
  have h2 : (1 : 𝒪[K.carrier]) = -1 := by simpa using congrArg Units.val h
  have h3 : ((2 : ℕ) : 𝒪[K.carrier]) = 0 := by push_cast; linear_combination h2
  simpa using Nat.cast_eq_zero.mp h3

/-- ★★★★**退化していない** —— `ker(Art) ⊓ I_K` は `I_K` そのものではない。

★これが無いと N1 は「`I_K = I_K`」という空虚な等式でありうる。
`Art(I_K) = 𝒪_K^×`(`ArtinFilteredDatum.map_absInertia`、仮定ゼロ)と
`𝒪_K^×` が自明でないことから出る。

★★同時に、`I_K / ‾⁅Γ_K,Γ_K⁆ ≅ 𝒪_K^×` という**本当の内容**がここに見えている。 -/
theorem topCommutator_ne_absInertia (K : PAdicLocalField p) (D : ArtinInertiaDatum K) :
    (commutator K.absGal).topologicalClosure ≠ absInertia K := by
  intro h
  haveI := nontrivial_units K
  obtain ⟨u, hu⟩ := exists_ne (1 : (𝒪[K.carrier])ˣ)
  have h1 : u ∈ Subgroup.map D.art (absInertia K) :=
    D.toArtinFilteredDatum.map_absInertia ▸ Subgroup.mem_top u
  obtain ⟨g, hg, rfl⟩ := h1
  have hmem : g ∈ D.art.ker ⊓ absInertia K := by
    rw [D.ker_inf_absInertia, h]; exact hg
  exact hu hmem.1

def topCommutator_ne_absInertia.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

section Compare

variable (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    [Normal K.carrier (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf)]

/-- ★★★**同じ `K` の上の 2 つの Lubin-Tate データは、惰性群に制限した核が一致する。**

★★これは `ReciprocityAlphaTransport::artinKerTransport_refl_iff_lubinTateClosure_eq` が
「制限しない版では `K_π = K_{π′}` と同値」と言った、まさにその状況である。
`K_π = K_{π′}` は偽(反例 `K = ℚ_p`、`π = p`、`π′ = -p`)だが、
★**惰性群に制限すると無条件に一致する。** -/
theorem ker_inf_absInertia_eq_of_lubinTate
    {π' : 𝒪[K.carrier]} (hπ'max : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π'})
    (hπ'ne0 : π' ≠ 0)
    (f' : PowerSeries 𝒪[K.carrier]) (hf'0 : PowerSeries.coeff 0 f' = 0)
    (hf'1 : PowerSeries.coeff 1 f' = π')
    (hf' : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f' = PowerSeries.X ^ (pp ^ ff))
    [Normal K.carrier (lubinTateClosure K hq hπ'max hπ'ne0 f' hf'0 hf'1 hf')] :
    (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf).ker ⊓ absInertia K
      = (reciprocityUnits K hq hπ'max hπ'ne0 f' hf'0 hf'1 hf').ker ⊓ absInertia K := by
  rw [ker_reciprocityUnits_inf_absInertia_eq_topCommutator K hq hπmax hπne0 f hf0 hf1 hf,
    ker_reciprocityUnits_inf_absInertia_eq_topCommutator K hq hπ'max hπ'ne0 f' hf'0 hf'1 hf']

def ker_inf_absInertia_eq_of_lubinTate.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Theorem 6.15", sectionId := "thm-6-15" }

end Compare

#print axioms restrictEquiv
#print axioms surjective_comp_subtype
#print axioms map_restrictEquiv_ker
#print axioms transportOfKerInfEq
#print axioms transportOfKerInfEq_apply
#print axioms map_map_transportOfKerInfEq
#print axioms fixingSubgroup_fixedField_of_isClosed
#print axioms fixingSubgroup_abelianClosure
#print axioms fixingSubgroup_abelianClosure_le_absInertia
#print axioms ker_reciprocityUnits_inf_absInertia
#print axioms ker_reciprocityUnits_inf_absInertia_eq_topCommutator
#print axioms ker_reciprocityUnits_inf_absInertia_inf
#print axioms nonempty_artinInertiaDatum
#print axioms artinKerInertiaTransport
#print axioms map_absInertia_filteredIso
#print axioms unitsTransportInertia
#print axioms unitsTransportInertia_art
#print axioms exists_unitsEquiv_inertia
#print axioms unitsTransportInertia_unique
#print axioms map_principalUnits_unitsTransportInertia
#print axioms unitsQuotientTransportInertia
#print axioms residueUnitsTransportInertia
#print axioms integersTransportInertia
#print axioms nonempty_integers_addEquiv_of_filteredIso
#print axioms nontrivial_units
#print axioms topCommutator_ne_absInertia
#print axioms ker_inf_absInertia_eq_of_lubinTate

end ABC3.Found.PGC
