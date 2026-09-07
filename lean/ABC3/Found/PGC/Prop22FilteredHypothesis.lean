import ABC3.Found.PGC.Prop22IntegersSuffice
import ABC3.Found.PGC.RamificationFiltrationZero
import ABC3.Found.PGC.RamificationNaturality
import ABC3.Skeleton.PGC.Section3Defs

/-!
# [pGC] Proposition 2.2 —— **分岐フィルトレーションを配線した**形

## ★★1. なぜこのファイルが要るか(D33)

原文 Proposition 2.2 は「`Γ_K` **と** `Γ_K^v`(`v > 0`)が与えられたとき」と言う。
ところが `Skeleton/PGC/Section2Defs.lean::RecoverableAsAddModule` が量化する `α` は
**位相群の同型だけ**で、`Skeleton/PGC/Section2.lean::prop_2_2` の
`(_RF : RamificationFiltration p)` は**先頭 `_` の未使用引数**である。

⇒ `Found/PGC/AbsClosureModules.lean::IntKbarRecoverable` は
**原典より強い可能性がある**(Noether の定理より、野性分岐のある `L/K` では
`𝒪_L` は `𝒪_K[Gal(L/K)]`-自由ではない ⇒ 正規底経由の Prop 2.1 の議論は
`𝒪` の側にはそのままでは効かない)。★偽だとは示されていない。

本ファイルは ★**濾過つきの変種を足す**。`RecoverableAsAddModule` 自体は
**変えない**(Prop 2.1 は原典どおり「`Γ_K` だけ」から復元する主張であり、
`Found/PGC/CountableGenerators.lean::prop_2_1` で無条件に閉じている)。
量化の形は `Found/PGC/Section3RealParameters.lean::Cor31Pinned` / `Cor33Pinned` と
**同じ**にした —— `α : FilteredGroup.Iso (filteredGroupOf RF K) (filteredGroupOf RF K')`、
`RF` は実物 `ABC3.Found.PGC.ramificationFiltration p`
(`Found/PGC/UnramifiedBaseChangeInvariance.lean`、仮定ゼロ)に固定する。
★自由な `RF` を許すと退化した `Gv ≡ ⊤` が入ってきて意味が消える(D32)。

## ★★2. 抽象核 —— `PAdicLocalField` も `p` も 1 語も出ない

§1 は **`FilteredGroup` だけ**の話である。分岐・付値・Galois・ノルム・`ℚ_p` は
1 語も出てこない。中身は 2 行で言える:

* `FilteredGroup.Iso` は **群体**をなす(`filteredIsoRefl` / `filteredIsoSymm` /
  `filteredIsoTrans`)。
* **内部自己同型は常に濾過を保つ** —— `Gv v` が**正規**だからである
  (`filteredIsoInner`、中身は mathlib `Subgroup.Normal.map_conj_eq` 1 本)。

★2 番目が非空虚性を無条件に与える: 濾過を課しても `Iso` は空にならない。

## ★★3. 直前の波との接続(名指し)

`Found/PGC/Prop22IntegersSuffice.lean` の実装者は
「本ファイルの主結果はすべて **α ごと**の形なので、後で α にフィルトレーション
両立性を課しても**そのまま使える**」と書いた。★**確かめた。そのまま繋がった。**
具体的には次の 7 本が**書き換えなしで**濾過つきの形に持ち上がる:

| 直前の波の宣言 | 濾過つきの持ち上げ先 |
|---|---|
| `compKbar_transport_of_intKbarTransport` | `prop_2_2_filtered_of_intKbarFiltered` |
| `closure_transport_of_intKbarTransport` | `recoverableAsAddModuleFiltered_closure_of_intKbarFiltered` |
| `intKbarTransport_refl` | `intKbarTransportFiltered_refl` |
| `intKbarTransport_symm` | `intKbarTransportFiltered_symm` |
| `intKbarTransport_trans` | `intKbarTransportFiltered_trans` |
| `intKbarTransport_inner` | `intKbarTransportFiltered_inner` |
| `intKbarTransport_galContinuousMulEquiv` | `intKbarTransportFiltered_algEquiv` |

★接続に要ったのは `filteredIsoEquiv`(1 層の射影、`rfl`)だけである。
`filteredIsoEquiv_refl` / `_symm` / `_trans` / `_inner` は **4 本とも `rfl`**。

## ★★★4. 退化の判定 —— 「今度は逆に自明になっていないか」

★★**結論: 自明にはなっていない。ただし `v ≤ 0` の部分は完全に無内容である。**

1. ★**非空虚**: `FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K)` は
   恒等と**すべての内部自己同型**を含む(`filteredIsoRefl` / `filteredIsoInner`、無条件)。
   ⇒ `IntKbarRecoverableFiltered` は「仮説が空虚だから真」という形ではない。
2. ★**強すぎない**: `IntKbarRecoverable ⟹ IntKbarRecoverableFiltered`
   (`intKbarRecoverableFiltered_of_intKbarRecoverable`)。★濾過を課すのは
   仮説を**弱める**向きであって、新しい要求を足してはいない。
3. ★★**`v ≤ 0` の条件は無内容**(`map_ramificationFiltration_of_nonpos`)——
   `Γ_K^v = I_K`(`v ≤ 0`、`ramificationFiltration_Gv_of_nonpos`)と
   Corollary 1.3(`inertia_recoverable_real`、無条件)から、
   ★**どんな連続同型 `α` でも `α(Γ_K^v) = Γ_{K'}^v` が `v ≤ 0` で自動的に成り立つ**。
   ⇒ 濾過を配線して増えた制約は **`v > 0` の部分だけ**である。
   ★これは原典が `v > 0` しかデータに入れていないことと**ぴったり一致する**
   (`map_ramificationFiltration_iff_pos` / `filteredIsoOfPos` で型として言う)。
   ★すなわち ★**我々の `FilteredGroup.Iso`(全実数 `v` で条件を課す)は
   原典のデータ(`v > 0` だけ)より強くない** —— 逸脱ではない。
4. ★**退化フィルトレーションではない**: `Γ_K^0 = I_K ≠ ⊤`
   (`ramificationFiltration_Gv_zero_ne_top`、無条件)。
   `Check/PGC/Theorem42NaiveGC.lean::topIsoOfContinuousMulEquiv` が示すとおり、
   `Gv ≡ ⊤` を入れると `map_Gv` は `Subgroup.map α ⊤ = ⊤` に潰れて
   **どんな `α` も濾過つき同型になる**。★実物の `ramificationFiltration p` では
   そうならない。
5. ★★**未解決として残ること(測っていないのではなく、測れていない)**:
   `v > 0` の条件が実際に `α` を切り落とすか——すなわち
   ★**濾過を保たない連続同型 `Γ_K ≃ₜ* Γ_{K'}` が存在するか**——は**分からない**。
   もし全ての連続同型が上付き分岐濾過を保つなら、`IntKbarRecoverableFiltered` は
   `IntKbarRecoverable` と**同値**になり、D33 の修理は空振りになる。
   ★これは「上付き分岐濾過が `Γ_K` から群論的に復元できるか」という問いそのもので、
   原典が `Out_Filt` を `Out` と区別している以上、原典も**開いたままにしている**。
   ★★**したがって D33 の修理は「効くかもしれないが、効くと示せてはいない」。**

## ★★5. 何が示せて、何が示せなかったか(生成元ごと)

群体なので生成元の `α` だけ見ればよい。★**濾過を課しても、`IntKbarTransport` が
言える `α` の集合は直前の波から 1 つも増えていない**:

| 生成元 | 濾過つきの同型か | `IntKbarTransport` |
|---|---|---|
| 恒等 | ○(無条件) | ○(無条件) |
| 内部自己同型 `g ↦ cgc⁻¹` | ○(無条件、正規性だけ) | ○(無条件) |
| 体の同型 `β : K ≃ₐ[ℚ_p] K'` から来る `α` | △(`IsNaturalFiltration` を仮定) | ○(無条件) |
| それ以外(Jarden–Ritter 型) | 不明 | ★**未解決**(=残る穴) |

★**残る穴は移動していない**: 「体の同型から来ない `α`」に対する
`IntKbarTransport` が依然として唯一の穴である。濾過はその穴を**塞いでいない**が、
穴の**入口を狭めうる**(上の 5 のとおり、狭まると示せてはいない)。

★次のノードの候補を名指しする:
`Found/PGC/PadicLogIntegers.lean::smul_padicLog_image_ramificationFiltration_eq_integers`
(`p^{-r}·log(Art(Γ^{re}_K)) = 𝒪_K`)は「`Γ^v` から `𝒪_K` を作る」を仮定ゼロで持っている。
これを `α` に沿って移送するには `reciprocityUnits` と `padicLog` の `α`-同変性が要る
——本シフトでは着手していない(Lubin–Tate のデータ `π`・`f` の選択に依存する形なので、
まず選択非依存性を要する)。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. `Skeleton/**` は 1 行も書き換えていない。`RecoverableAsAddModule` も
   `IntKbarRecoverable` も**消していない**。本ファイルは**変種を足しただけ**である。
2. `FilteredGroup.Iso` は全実数 `v` で `map_Gv` を課すが、原典が与えるデータは
   `v > 0` だけである。この差が**無い**ことを `map_ramificationFiltration_iff_pos` で
   証明した(上の 4-3)。★したがってこの点は逸脱では**ない**。
3. `filteredIsoOfAlgEquiv` / `intKbarTransportFiltered_algEquiv` は
   `IsNaturalFiltration (ramificationFiltration p)`(**未証明**)を明示的な仮説に取る。
   `Found/PGC/RamificationNaturality.lean` と同じ扱いである。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC

/-! ## §1 ★★抽象核 —— `FilteredGroup` だけの話

★分岐・付値・Galois・ノルム・`ℚ_p`・`PAdicLocalField` は 1 語も出てこない。
出てくるのは「位相群」「その正規閉部分群の反単調族」だけである。 -/

section AbstractCore

/-- **恒等**は濾過を保つ。 -/
def filteredIsoRefl (A : FilteredGroup) : FilteredGroup.Iso A A where
  equiv := ContinuousMulEquiv.refl A.G
  map_Gv v := by
    ext g
    exact ⟨by rintro ⟨h, hh, rfl⟩; exact hh, fun hg => ⟨g, hg, rfl⟩⟩

/-- **合成**は濾過を保つ。 -/
def filteredIsoTrans {A B C : FilteredGroup} (f : FilteredGroup.Iso A B)
    (g : FilteredGroup.Iso B C) : FilteredGroup.Iso A C where
  equiv := f.equiv.trans g.equiv
  map_Gv v := by
    rw [← g.map_Gv v, ← f.map_Gv v, Subgroup.map_map]
    rfl

/-- **逆**は濾過を保つ。`map_Gv` が「像がちょうど一致する」形なので出る
(包含だけだと出ない)。 -/
def filteredIsoSymm {A B : FilteredGroup} (f : FilteredGroup.Iso A B) :
    FilteredGroup.Iso B A where
  equiv := f.equiv.symm
  map_Gv v := by
    have hid : ((f.equiv.symm.toMulEquiv : B.G ≃* A.G) : B.G →* A.G).comp
        ((f.equiv.toMulEquiv : A.G ≃* B.G) : A.G →* B.G) = MonoidHom.id A.G :=
      MonoidHom.ext (fun x => f.equiv.symm_apply_apply x)
    rw [← f.map_Gv v, Subgroup.map_map, hid, Subgroup.map_id]

/-- ★★**内部自己同型は常に濾過を保つ**。中身は「`Gv v` が正規である」ことだけ
(`FilteredGroup.isNormal` + mathlib `Subgroup.Normal.map_conj_eq`)。

★これが濾過つきの `Iso` の**非空虚性**を無条件に与える。 -/
noncomputable def filteredIsoInner (A : FilteredGroup) (c : A.G) : FilteredGroup.Iso A A where
  equiv :=
    { toFun := fun g => c * g * c⁻¹
      invFun := fun g => c⁻¹ * g * c
      left_inv := fun g => by group
      right_inv := fun g => by group
      map_mul' := fun a b => by group
      continuous_toFun := by fun_prop
      continuous_invFun := by fun_prop }
  map_Gv v := by
    have hn : (A.Gv v).Normal := A.isNormal v
    exact Subgroup.Normal.map_conj_eq (H := A.Gv v) c

@[simp] theorem filteredIsoInner_equiv_apply (A : FilteredGroup) (c g : A.G) :
    (filteredIsoInner A c).equiv g = c * g * c⁻¹ := rfl

/-- 濾過つき同型は空でない(どんな `FilteredGroup` でも)。 -/
theorem nonempty_filteredIso_self (A : FilteredGroup) :
    Nonempty (FilteredGroup.Iso A A) := ⟨filteredIsoRefl A⟩

end AbstractCore

variable {p : ℕ} [Fact p.Prime]

/-! ## §2 実物の分岐フィルトレーションに固定する

★`Cor31Pinned` / `Cor33Pinned`(`Found/PGC/Section3RealParameters.lean`)と
**同じ形**。`RF` は自由に量化せず、実物
`ramificationFiltration p`(`Found/PGC/UnramifiedBaseChangeInvariance.lean`、仮定ゼロ)
に固定する。 -/

/-- 実物の分岐フィルトレーションを載せた `Γ_K`。 -/
noncomputable def pgcFilteredGroup (K : PAdicLocalField p) : FilteredGroup :=
  filteredGroupOf (ramificationFiltration p) K

def pgcFilteredGroup.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Definition 2.3 (FilteredGroup.Iso)",
    sectionId := "def-2-3" }

/-- 濾過つき同型の台となる位相群同型。★1 層の射影なので `rfl` で足りる
(2 層をまたぐ `rfl` は kernel を止める —— `lean-idioms.md` #59)。 -/
def filteredIsoEquiv {K K' : PAdicLocalField p}
    (α : FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K')) :
    ContinuousMulEquiv K.absGal K'.absGal := α.equiv

@[simp] theorem filteredIsoEquiv_apply {K K' : PAdicLocalField p}
    (α : FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K')) (g : K.absGal) :
    filteredIsoEquiv α g = α.equiv g := rfl

/-! ## §3 濾過つきの回復可能性 -/

/-- ★★**濾過つきの `RecoverableAsAddModule`**。

`Skeleton/PGC/Section2Defs.lean::RecoverableAsAddModule` との違いは 1 点だけ ——
量化する `α` が**位相群の同型**ではなく**濾過つきの同型**であること。
★これが原文 Proposition 2.2 の「`Γ_K` **と** `Γ_K^v` が与えられたとき」に対応する。 -/
def RecoverableAsAddModuleFiltered (Obj : PAdicLocalField p → Type*)
    [∀ K, AddCommGroup (Obj K)] [∀ K, DistribMulAction K.absGal (Obj K)] : Prop :=
  ∀ {K K' : PAdicLocalField p}
    (α : FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K')),
    ∃ φ : Obj K ≃+ Obj K', ∀ (g : K.absGal) (x : Obj K),
      φ (g • x) = ((filteredIsoEquiv α).toMulEquiv g) • (φ x)

def RecoverableAsAddModuleFiltered.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- ★**濾過を課すのは仮説を弱める向き**。強くはしていない。 -/
theorem recoverableAsAddModuleFiltered_of_recoverableAsAddModule
    (Obj : PAdicLocalField p → Type*)
    [∀ K, AddCommGroup (Obj K)] [∀ K, DistribMulAction K.absGal (Obj K)]
    (h : RecoverableAsAddModule (p := p) Obj) : RecoverableAsAddModuleFiltered (p := p) Obj :=
  fun α => h (filteredIsoEquiv α)

/-- 濾過つきの `IntKbarRecoverable`。 -/
def IntKbarRecoverableFiltered (p : ℕ) [Fact p.Prime] : Prop :=
  RecoverableAsAddModuleFiltered (p := p) (fun K => IntKbar K)

def IntKbarRecoverableFiltered.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- 濾過つきの `CompKbarRecoverable`。 -/
def CompKbarRecoverableFiltered (p : ℕ) [Fact p.Prime] : Prop :=
  RecoverableAsAddModuleFiltered (p := p) (fun K => CompKbar K)

def CompKbarRecoverableFiltered.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- **α ごとの濾過つき `𝒪_{K̄}` の移送**。 -/
def IntKbarTransportFiltered {K K' : PAdicLocalField p}
    (α : FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K')) : Prop :=
  IntKbarTransport (filteredIsoEquiv α)

theorem intKbarRecoverableFiltered_iff :
    IntKbarRecoverableFiltered p ↔
      ∀ {K K' : PAdicLocalField p}
        (α : FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K')),
        IntKbarTransportFiltered α := Iff.rfl

theorem intKbarRecoverableFiltered_of_intKbarRecoverable (h : IntKbarRecoverable (p := p)) :
    IntKbarRecoverableFiltered p :=
  recoverableAsAddModuleFiltered_of_recoverableAsAddModule _ h

/-! ### 直前の波(`Prop22IntegersSuffice.lean`)がそのまま繋がる

★実装者の見立て「主結果はすべて **α ごと**の形なのでそのまま使える」は正しかった。
以下の 2 本は、直前の波の定理を**書き換えずに**適用しただけである。 -/

/-- ★★★**[pGC] Proposition 2.2(濾過を配線した形)—— 仮説は `𝒪_{K̄}` だけでよい**。 -/
theorem prop_2_2_filtered_of_intKbarFiltered (h : IntKbarRecoverableFiltered p) :
    IntKbarRecoverableFiltered p ∧ CompKbarRecoverableFiltered p :=
  ⟨h, fun α => compKbar_transport_of_intKbarTransport (α := filteredIsoEquiv α) (h α)⟩

def prop_2_2_filtered_of_intKbarFiltered.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- ★**Proposition 2.1(濾過つき)は `𝒪_{K̄}` の回復から従う**。 -/
theorem recoverableAsAddModuleFiltered_closure_of_intKbarFiltered
    (h : IntKbarRecoverableFiltered p) :
    RecoverableAsAddModuleFiltered (p := p) (fun K => K.closure) :=
  fun α => closure_transport_of_intKbarTransport (α := filteredIsoEquiv α) (h α)

/-! ## §4 群体性と非空虚性 —— 生成元の `α` だけ見ればよい

★`filteredIsoEquiv` を通した 4 本の等式は**すべて `rfl`** である。
だから直前の波の群体補題(`intKbarTransport_{refl,symm,trans}`)と
非空虚性(`intKbarTransport_inner` / `intKbarTransport_galContinuousMulEquiv`)が
そのまま持ち上がる。 -/

theorem filteredIsoEquiv_refl (K : PAdicLocalField p) :
    filteredIsoEquiv (filteredIsoRefl (pgcFilteredGroup K))
      = ContinuousMulEquiv.refl K.absGal := rfl

theorem filteredIsoEquiv_inner (K : PAdicLocalField p) (c : K.absGal) :
    filteredIsoEquiv (filteredIsoInner (pgcFilteredGroup K) c) = innerAbsGalEquiv K c := rfl

theorem filteredIsoEquiv_symm {K K' : PAdicLocalField p}
    (α : FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K')) :
    filteredIsoEquiv (filteredIsoSymm α) = (filteredIsoEquiv α).symm := rfl

theorem filteredIsoEquiv_trans {K K' K'' : PAdicLocalField p}
    (α : FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K'))
    (β : FilteredGroup.Iso (pgcFilteredGroup K') (pgcFilteredGroup K'')) :
    filteredIsoEquiv (filteredIsoTrans α β)
      = (filteredIsoEquiv α).trans (filteredIsoEquiv β) := rfl

theorem intKbarTransportFiltered_refl (K : PAdicLocalField p) :
    IntKbarTransportFiltered (filteredIsoRefl (pgcFilteredGroup K)) :=
  intKbarTransport_refl K

/-- ★**内部自己同型では無条件に成り立つ**(濾過つきでも同じ)。 -/
theorem intKbarTransportFiltered_inner (K : PAdicLocalField p) (c : K.absGal) :
    IntKbarTransportFiltered (filteredIsoInner (pgcFilteredGroup K) c) :=
  intKbarTransport_inner K c

theorem intKbarTransportFiltered_symm {K K' : PAdicLocalField p}
    {α : FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K')}
    (h : IntKbarTransportFiltered α) : IntKbarTransportFiltered (filteredIsoSymm α) :=
  intKbarTransport_symm h

theorem intKbarTransportFiltered_trans {K K' K'' : PAdicLocalField p}
    {α : FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K')}
    {β : FilteredGroup.Iso (pgcFilteredGroup K') (pgcFilteredGroup K'')}
    (hα : IntKbarTransportFiltered α) (hβ : IntKbarTransportFiltered β) :
    IntKbarTransportFiltered (filteredIsoTrans α β) :=
  intKbarTransport_trans hα hβ

/-- ★濾過つきの `Iso` は(自己同型では)**無条件に**空でない。 -/
theorem nonempty_filteredIso_pgc (K : PAdicLocalField p) :
    Nonempty (FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K)) :=
  nonempty_filteredIso_self _

/-- 自然性を仮定すると、体の同型からも濾過つき同型が作れる。

★`IsNaturalFiltration (ramificationFiltration p)` は**未証明**
(`Skeleton/PGC/Section4.lean` の注記と同じ扱い)。 -/
noncomputable def filteredIsoOfAlgEquiv (hnat : IsNaturalFiltration (ramificationFiltration p))
    {K K' : PAdicLocalField p} (β : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) :
    FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K') :=
  naturalFilteredIso (ramificationFiltration p) hnat β

/-- ★体の同型から来る `α` では、濾過つきでも**無条件に**成り立つ
(`hnat` は `α` を作るためだけに要る)。 -/
theorem intKbarTransportFiltered_algEquiv (hnat : IsNaturalFiltration (ramificationFiltration p))
    {K K' : PAdicLocalField p} (β : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) :
    IntKbarTransportFiltered (filteredIsoOfAlgEquiv hnat β) :=
  intKbarTransport_galContinuousMulEquiv β

/-! ## §5 ★★★退化の自己検査

「濾過を配線した形が、今度は逆に弱すぎ(自明)になっていないか」を型で調べる。 -/

/-- ★★★**`v ≤ 0` での `map_Gv` 条件は、どんな連続同型でも自動的に成り立つ**。

`Γ_K^v = I_K`(`v ≤ 0`、`ramificationFiltration_Gv_of_nonpos`)と
Corollary 1.3(`inertia_recoverable_real`、無条件)を合わせるだけ。

⇒ ★**濾過を配線して増えた制約は `v > 0` の部分だけ**である。 -/
theorem map_ramificationFiltration_of_nonpos {K K' : PAdicLocalField p}
    (α : ContinuousMulEquiv K.absGal K'.absGal) {v : ℝ} (hv : v ≤ 0) :
    Subgroup.map α.toMulEquiv ((ramificationFiltration p).Gv K v)
      = (ramificationFiltration p).Gv K' v := by
  rw [ramificationFiltration_Gv_of_nonpos K hv, ramificationFiltration_Gv_of_nonpos K' hv,
    absInertia_eq_inertia, absInertia_eq_inertia]
  exact inertia_recoverable_real α

def map_ramificationFiltration_of_nonpos.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Corollary 1.3", sectionId := "cor-1-3" }

/-- ★★**我々の `FilteredGroup.Iso` は原典のデータより強くない**。

`FilteredGroup.Iso` は全実数 `v` で `map_Gv` を課すが、原典が与えるのは
`Γ_K^v`(`v > 0`)だけである。★その差は**無い**——`v ≤ 0` の条件は自動だから。 -/
theorem map_ramificationFiltration_iff_pos {K K' : PAdicLocalField p}
    (α : ContinuousMulEquiv K.absGal K'.absGal) :
    (∀ v : ℝ, Subgroup.map α.toMulEquiv ((ramificationFiltration p).Gv K v)
        = (ramificationFiltration p).Gv K' v)
      ↔ (∀ v : ℝ, 0 < v → Subgroup.map α.toMulEquiv ((ramificationFiltration p).Gv K v)
        = (ramificationFiltration p).Gv K' v) := by
  refine ⟨fun h v _ => h v, fun h v => ?_⟩
  rcases lt_or_ge 0 v with hv | hv
  · exact h v hv
  · exact map_ramificationFiltration_of_nonpos α hv

def map_ramificationFiltration_iff_pos.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- ★★**原典のデータ(`v > 0` だけ)から濾過つき同型を作る**。
`map_ramificationFiltration_iff_pos` の構成的な形。 -/
noncomputable def filteredIsoOfPos {K K' : PAdicLocalField p}
    (α : ContinuousMulEquiv K.absGal K'.absGal)
    (h : ∀ v : ℝ, 0 < v → Subgroup.map α.toMulEquiv ((ramificationFiltration p).Gv K v)
      = (ramificationFiltration p).Gv K' v) :
    FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K') where
  equiv := α
  map_Gv v := (map_ramificationFiltration_iff_pos α).mpr h v

/-- ★★**自然性も `v > 0` だけ確かめれば十分**。
`Found/PGC/RamificationNaturality.lean::IsNaturalFiltration` を実物に固定した形。
★次のノードの証明義務が(実質)半分になる。 -/
theorem isNaturalFiltration_ramificationFiltration_iff_pos :
    IsNaturalFiltration (ramificationFiltration p) ↔
      ∀ {K K' : PAdicLocalField p} (β : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (v : ℝ), 0 < v →
        Subgroup.map (galContinuousMulEquiv β).toMulEquiv ((ramificationFiltration p).Gv K v)
          = (ramificationFiltration p).Gv K' v := by
  constructor
  · intro h _ _ β v _
    exact h β v
  · intro h _ _ β v
    exact (map_ramificationFiltration_iff_pos (galContinuousMulEquiv β)).mpr
      (fun w hw => h β w hw) v

/-- ★★**実物のフィルトレーションは退化していない**: `Γ_K^0 = I_K ≠ ⊤`。

`Check/PGC/Theorem42NaiveGC.lean::topIsoOfContinuousMulEquiv` が示すとおり、
`Gv ≡ ⊤` を入れると `map_Gv` は `Subgroup.map α ⊤ = ⊤` に潰れて
**どんな `α` も濾過つき同型になる**。★実物ではそうならない。 -/
theorem ramificationFiltration_Gv_zero_ne_top (K : PAdicLocalField p) :
    (ramificationFiltration p).Gv K 0 ≠ ⊤ := by
  rw [ramificationFiltration_Gv_zero K]
  exact absInertia_ne_top K

/-- ★同じことを `pgcFilteredGroup` の言葉で。 -/
theorem pgcFilteredGroup_Gv_zero_ne_top (K : PAdicLocalField p) :
    (pgcFilteredGroup K).Gv 0 ≠ ⊤ := ramificationFiltration_Gv_zero_ne_top K

#print axioms filteredIsoRefl
#print axioms filteredIsoSymm
#print axioms filteredIsoTrans
#print axioms filteredIsoInner
#print axioms prop_2_2_filtered_of_intKbarFiltered
#print axioms recoverableAsAddModuleFiltered_closure_of_intKbarFiltered
#print axioms intKbarTransportFiltered_inner
#print axioms map_ramificationFiltration_of_nonpos
#print axioms map_ramificationFiltration_iff_pos
#print axioms filteredIsoOfPos
#print axioms isNaturalFiltration_ramificationFiltration_iff_pos
#print axioms ramificationFiltration_Gv_zero_ne_top

end ABC3.Found.PGC
