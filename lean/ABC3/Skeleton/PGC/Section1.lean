import ABC3.Skeleton.PGC.Setup
import ABC3.Interface.PGC.LocalFieldData

/-!
# [pGC] §1 — 命題

設定・記号・§1冒頭の暗黙の定義は `ABC3/Skeleton/PGC/Setup.lean`。
未構築の基礎は `ABC3/Interface/PGC/LocalFieldData.lean`。
-/

namespace ABC3.Skeleton.PGC

open ABC3.Meta ABC3.Interface.PGC

variable {p : ℕ} [Fact p.Prime]

/-! ## Proposition 1.1 -/

/-- 円分指標 `χ_K : Γ_K → ℤ_[p]ˣ` を「K に付随する対象」として束ねたもの。

`χ_K` 自体は mathlib の `cyclotomicCharacter` から得る——**これは K(と K̄)を使って定義される**。
Proposition 1.1 が主張するのは、こうして作った対象が実は Γ_K だけで決まる、ということ。
移送は「α で引き戻す」= `f ∘ α.symm`。 -/
noncomputable def cyclotomicCharacterObject : AssociatedObject p where
  Obj := fun K => K.absGal → ℤ_[p]ˣ
  obj := fun K g => cyclotomicCharacter K.closure p g.toRingEquiv
  transport := fun α f g' => f (α.toMulEquiv.symm g')

def cyclotomicCharacterObject.src : Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- **[pGC] Proposition 1.1**

原文 (pGC p.3):
> The cyclotomic character χ : Γ_K → Z[bb]_p^× can be recovered entirely
> group-theoretically from Γ_K.

## 未解決(1_Structured の `<p class="open">` と対応)

原典は Proposition 1.1 に独立した証明を付けておらず、直前の地の文(双対性の復習)が論拠になる。
しかしその地の文が確立するのは「Γ_K-加群 **Z**_p(1) の**同型類**が回復できる」ことであり、
本命題が主張するのは「**指標** χ が回復できる」こと。**この間の一段は原典にない。**
証明を付ける段階で、
  ① 我々のモデル化の誤り(§1冒頭の定義から直ちに従う)か
  ② 必要な数学が未構築(局所 Tate 双対性——mathlib 不在)か
  ③ 原典側の飛躍か
を切り分けること。既定は①。

## 依拠する境界外の結果(Track B のキュー候補)

- [2] Bloch–Kato, Grothendieck Festschrift I (1990), Proposition 3.8 — 局所 Tate 双対性。
  mathlib 不在(実測 v4.31.0-rc2)。

## 反証可能性(2026-08-14 監査)

Corollary 3.12 が現在の `Interface` の下で偽だったこと(`Check/IUTchIII/Cor312Degenerate.lean`)を
受けて横断監査した。**この定理は反証できなかった。** 本定理は `Interface` を取らないうえ、
共役による経路は `Check.PGC.cyclotomicCharacter_conj_invariant` で
**閉じていることが証明済み**(円分指標は可換群への準同型なので内部自己同型で動かない)。
探した範囲は `Check/PGC/RefutationAttempts.lean`。 -/
theorem cyclotomicCharacter_recoverable :
    (cyclotomicCharacterObject (p := p)).RecoverableFromAbsGal := by
  sorry

def cyclotomicCharacter_recoverable.src : Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- 原文の証明文から抽出した、証明が要求するもの(G6)。

★これは **下界** である——原文が挙げていない依存は写らない。 -/
def cyclotomicCharacter_recoverable.needs : List ProofObligation :=
  [ .citation "[2] Bloch-Kato, Grothendieck Festschrift I (1990)"
      "Proposition 3.8(局所 Tate 双対性)"
      (.absent <|
        "mathlib v4.31.0-rc2 全体を再測定(2026-08-14): " ++
        "grep -rlE 'localTateDuality|local Tate|TateDuality|tateDuality' → 0 件、" ++
        "grep -rlE '\\bTate\\b' → 4 件だが全て別物" ++
        "(EllipticCurve/Weierstrass・RingTheory/Adjoin/Tower・Perfectoid/BDeRham・Perfectoid/FontaineTheta)、" ++
        "grep -rlE 'tateCohomology|TateCohomology|GaloisCohomology' → 0 件、" ++
        "grep -rli 'duality' Mathlib/NumberTheory Mathlib/RepresentationTheory → " ++
        "DirichletCharacter/Orthogonality・MulChar/Duality・Cyclotomic/Galois・Tannaka のみ。" ++
        "公開プロジェクト(ClassFieldTheory・LocalClassFieldTheory)は 2026-08-14 に測定して 0 ヒット" ++
        "(ResearchPaper/lean-ecosystem.json。両リポジトリは本リポジトリに clone されていないので今日は再現していない)。" ++
        "★★2026-09-04 再測定(mathlib が更新されていた——上の 2026-08-14 の" ++
        "'tateCohomology|TateCohomology|GaloisCohomology → 0件' は今は誤り、古い版の話):" ++
        "`RepresentationTheory/Homological/GroupCohomology/`(groupCohomology・H0/H1/H2、" ++
        "有限群の非連続コホモロジー)と `RepresentationTheory/Homological/TateCohomology/` が" ++
        "存在する。ただし `groupCohomology` は `[Group G]` のみを要求する**離散**群論で、" ++
        "副有限 Γ_K の連続コサイクルではない——有限次拡大 L/K の Gal(L/K)(有限群)には" ++
        "直接使える可能性があるが、(a) H²(Gal(L/K),Lˣ) と Brauer 群 Br(L/K) を繋ぐ" ++
        "crossed product 構成、(b) 不変写像 Br(K)≅ℚ/ℤ の具体形、(c) Γ_K(射影極限)への" ++
        "extension のいずれも mathlib に無い(`BrauerGroup`(`Algebra/BrauerGroup/`)は" ++
        "CSA 経由の抽象論のみで、局所体の Brauer 群の計算は無い)。**抽象論の部品は増えたが、" ++
        "局所体特有の算術内容(不変写像の具体形)は依然として不在**——次にここへ戻るなら" ++
        "(a)(b)(c) のどれか1つを先に測ること") 2,
    .derivation "H^2(K,M) ≅ Z/p^nZ ⟺ M ≅ Z/p^nZ(1) を双対性から導く段" 3,
    .implicitStep
      "★★★★★★★★★★2026-09-06 の再測定(第 1038)。★★★**2026-09-06 後刻の訂正: この欄の旧記述「H2 経由は捨ててよい。H1 経由の在庫がある」は**楽観的すぎた**。`H1(K, mu_n) = K^x/(K^x)^n` は確かに在庫だが、★**Kummer 理論は `mu_n` を係数として最初から持っている**ので、**円分子の構成には使えない**(構成したいものを入力に要求している)。H1 が H2 を代替できるのは局所 Euler 標数公式による計数経路だけで、その上界の標準証明は双対性を通る。★★**円分子は必ず局所相互律か局所双対性のどちらかを通る**。安い迴回路は 3 つとも潰した(Kummer 単独 / 不分岐指標 / H1 の計数)。推奨は**経路 Λ(LCFT 円分子)**で、`tors_{p^n}(Gal_L^{ab}) = mu_{p^n}` を使う道。見積は 11 ノード・中央値 8,700 行。**`Found/PGC/KummerDuality.lean` の `KummerDual.range_kummerMap` と `KummerDual.ker_kummerMap` を合わせると **F^x/(F^x)^n ≅ Hom_cont(Γ_F, μ_n)**、すなわち **H1(K, μ_n) ≅ K^x/(K^x)^n** であり、Bloch-Kato Prop 3.8 の**次数 1 の半分**である(mathlib 不在のまま自前に証明済み、sorry 無し)。★さらに `Found/PGC/LubinTateReciprocityIsomorphism.lean:404` の **`galoisReciprocityEquiv : Gal(K⟮ x⟯/K) ≃* (𝒪_K ⧸ π^n)^x`** がある——**局所相互律は mathlib に全く不在**(re:artinMap|ArtinMap|reciprocityIso|localReciprocity → 0)なのに、この木は Lubin-Tate 経由で自前に持っている。★本ファイルは `Skeleton.PGC.Setup` と `Interface.PGC.LocalFieldData` しか import しておらず、**同じ木の中の完成品を見ていなかった**(「`Found` に在るのに `Skeleton` が参照していない」の **5 件目**)" 3,
    .citation "[mathlib]" "連続コホモロジー((c) Γ_K への extension)"
      (.inMathlib "continuousCohomology") 3,
    .implicitStep
      "★(c) の使い方の注意(2026-09-06 実測)。`continuousCohomology (R G) (n) : Action (TopModuleCat R) G ⥤ TopModuleCat R` (`Mathlib.Algebra.Category.ContinuousCohomology.Basic`、216 行)。★対象が `Rep R G` ではなく **`Action (TopModuleCat R) G`** なので橋が 1 本要る。★**長完全列が無い**(TODO)し、`groupCohomology` との一致も TODO。つまり「定義はある、道具はまだ無い」。★★**2026-09-04 の記録が (c) を「不在」としていたのは誤り**である" 3,
    .citation "[mathlib]" "(a) crossed product / (b) 不変写像 Br(K) の具体形 / 局所 Tate 双対性"
      (.absent
        ("2026-09-06 再測(★索引の欠陥を直し、非 ASCII を含む名前 24,659 件が載った状態での測定): "
         ++ "re:crossed → 0 / re:crossedProduct|CrossedProduct → 0 / re:RelativeBrauer → 0 / "
         ++ "re:invariantMap|BrauerInvariant|invBrauer|hasseInvariant → 0 / "
         ++ "re:normResidueSymbol|normResidue → 0 / "
         ++ "re:artinMap|ArtinMap|reciprocityIso|localReciprocity → 0 / "
         ++ "re:localTateDuality|TateDuality|tateDuality|PoitouTate|Poitou → 0 / re:Herbrand → 0。 "
         ++ "★re:Brauer → 9 件だが**全 9 件が Algebra/BrauerGroup/Defs.lean(全 100 行)**で、"
         ++ "CSA の商集合があるだけ。**BrauerGroup は Type であって群ですらない**"
         ++ "(ファイル冗頭の TODO 1 が「アーベル群であることを証明する」)。 "
         ++ "★左辺(H2 側)は書ける——`Rep.ofAlgebraAutOnUnits` と `groupCohomology.H2` がある。"
         ++ "**繋ぐ側だけが無い**。 "
         ++ "★`tateCohomology` は 33 件あるが **[Fintype G] を要求する**ので Γ_K には使えない"
         ++ "(有限次 Gal(L/K) 専用。ファイル冗頭に 2025 ClassFieldTheory workshop の成果と明記)。 "
         ++ "★Tate コホモロジーの慣習的な綴り re:Ĥ も **0 件**"
         ++ "(名前欄・statement 欄の両方、Ĥ|Ħ|Ḣ の全異体で 0)——"
         ++ "★索引の切れではなく**本当に不在**。mathlib は ASCII で tateCohomology と綴る")) 2,
    .implicitStep "Z_p(1) の同型類が回復できることから、指標 χ が回復できることへの一段" 3,
    .implicitStep
      "★★★★★★2026-09-06(第 1038): **docstring の ①②③ の切り分けが決着した**。『`Z_p(1)` の同型類 ⇒ 指標 χ』の一段は **`Found/PGC/CyclotomicRecovery.lean::padicUnits_eq_of_smul_equivariant` で閉じた**(階数 1 の `Z_p`-表現は同型類が指標を決める。`Z_p`-線形同型は `φ 1` 倍しかなく、整域で約せる。sorry 無し)。★★したがって **③(原典の飛躍)ではない**。残るのは ② だけである" 3,
    .implicitStep
      "★★★★2026-09-06(第 1038): **α が体の同型から来る場合は無条件に成立する**(`cyclotomicCharacterObject_transport_galContinuousMulEquiv`、sorry 無し)。鍵は `cyclotomicCharacter_conj`(**円分指標の自然性**、mathlib に不在)である。★この形式化が展開すると厳密に `∀ K K' (α : Γ_K ≃ₜ* Γ_K') (g), χ_{K'}(α g) = χ_K(g)` になることを Lean で確定させた(`cyclotomicCharacterObject_recoverable_iff`)。★★これは §1 冒頭の原典自身の定義そのままで、Prop 1.2 の `∀ RD` のような**自由なデータ引数が無い**ので、余計な強さが混入する余地がない(★D13 で Prop 1.2 に起きた「修復が強すぎる主張を作った」はここでは起きていない)" 3,
    .implicitStep
      "★★★★★★★★**残る壁は 1 つだけ**(2026-09-06)。Jarden-Ritter 型の α(体の同型から来ない α)に対しては、**円分子 `Λ(Γ_K) = Hom(H2(Γ_K, Z/p^n), Z/p^n) ≅ μ_{p^n}`** の群論的構成が要る。接続点は `cyclotomicCharacterObject_transport_of_moduleEquiv` に置いてあり、その仮説がそのまま次のノードの statement になる。★★**前回の段取りは誤りだった** —— `galoisReciprocityEquiv` を `K = Q_p` に特化して `Gal(K(μ_{p^n})/K) ≃ (Z/p^n)^x` を得ても Prop 1.1 には**届かない**。`RecoverableFromAbsGal` は**任意の抽象位相群同型 α** についての主張なので、各 K ごとの `Gal(K(μ_{p^n})/K)` の構造を知っても α の振る舞いは制約されない" 3,
    .citation "[この木]" "`iteratedLubinTatePsiTorsionPoints` と `μ_{p^n}` を同定する橋"
      (.absent
        ("2026-09-06 実測: `.cache/decl-index.txt` で "
         ++ "re:iteratedLubinTatePsiTorsionPoints → 14 件を**全走査**したが、"
         ++ "re:rootsOfUnity|IsPrimitiveRoot に触れるものは **0 件**。"
         ++ "`Found/PGC/QpRootsOfUnity.lean` は `l ≠ p` の素な乗根の話で p 冪乗根を扱っていない。"
         ++ "★Prop 1.1 には直結しないが、`galoisReciprocityEquiv` を円分側で使うための"
         ++ "独立に有用な部品である")) 3 ]


/-! ## Proposition 1.2 -/

/-- Proposition 1.2 が回復されると主張する2つの量を「K に付随する対象」として束ねたもの:
剰余体の元の個数 q と、絶対次数 [K : ℚ_p]。

いずれも**数**なので、移送は恒等——原典の「α によって対応物へ移る」は、数については
「等しい」を意味する。

`q` は `Interface` の `ResidueCardinality` から取る(まだ構築できていないため)。
`[K : ℚ_p]` は mathlib の `Module.finrank` で書ける。 -/
noncomputable def residueCardAndDegreeObject (RD : ResidueCardinality p) :
    AssociatedObject p where
  Obj := fun _ => ℕ × ℕ
  obj := fun K => (RD.card K, Module.finrank ℚ_[p] K.carrier)
  transport := fun _ x => x

def residueCardAndDegreeObject.src : Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.2", sectionId := "prop-1-2" }

/-- **[pGC] Proposition 1.2**

原文 (pGC p.3):
> The number q of elements in the residue field of O[scr]_K, and well as the absolute
> degree [K : Q[bb]_p] of K, can be recovered entirely group-theoretically from Γ_K.

原文の "and well as" は "as well as" の誤植と読めるが、逐語欄では原文どおりに保つ。

## 条件付き形式化

`RD : ResidueCardinality p` を仮説として取る——q を与える構成そのものが未構築だから
(`Interface/PGC/LocalFieldData.lean`)。`RD` が本物になれば、この定理は無条件になる。

## 未解決

原典の論拠は局所類体論の同型 Γ_K^ab ≅ (K^×)^∧ と、そこからの計数
(捩れの prime-to-p 部分が q−1 個、pro-p 商の階数が [K:ℚ_p]+1)。
後者は p進対数を使う——**mathlib にも公開プロジェクトにも無い**(実測)。

## 反証可能性(2026-08-14 監査)

**反証できなかった。** `RD : ResidueCardinality p` が自由なデータであることは
**この定理を反証する役に立たない**——K′ = K の経路は `transport` が恒等なので
自明に閉じ(`Check.PGC.residueCardAndDegree_self_route_closed`)、
残るのは「異なる2体の間の連続同型 α」の構成だけ。しかも
`Check.PGC.refutation_reduces_to_alpha` が示すとおり、その α が1つでも作れれば
`RD` に関係なく次数成分だけで落ちる——すなわち**反証するには本命題自身を偽にするしかない**。
探した範囲は `Check/PGC/RefutationAttempts.lean`。

## ★2026-09-05: 上の監査は誤りだった——修理の記録(7 例目)

`Check/PGC/Prop12Degenerate.lean`(第 1012)が反証を作った。
上の監査は「異なる2体の間の α が要る」までは正しかったが、
**その α は体の同型から作れる**ことを見落としていた——台の型を `ℚ_[p]` のままにして
体構造だけを `x ↦ -x` に沿って移送すると、ℚ_p-代数として同型なのに
`PAdicLocalField p` の項としては異なる2体ができ、
`Found/PGC/GaloisTransferContinuous.lean::galContinuousMulEquiv` が α を与える。
次数成分は両方 1 なので効かず、落ちるのは `card` 成分だった。

落とした条件は**同型不変性**——「落とした条件は主張を偽にするか自明にする」の 7 例目。
修理として `Interface/PGC/LocalFieldData.lean` の `ResidueCardinality` に
`card_congr`(ℚ_p-代数同型なら同じ q)を足した。実物がそれを満たすことは
`Found/PGC/ResidueCardinality.lean::residueCard_congr` で証明済み(`sorry` 無し)——
ℚ_p-代数同型はスペクトルノルムを保つので整数環の環同型を誘導し、剰余体の濃度が一致する。
修理後は反例が作れないことも `Prop12Degenerate.lean::no_residueCardinality_with_badRD_card`
で確認した。★これで本定理は原典本来の内容(α から体の同型を作る)に戻っており、
**易しくはなっていない**。 -/
theorem residueCard_and_degree_recoverable (RD : ResidueCardinality p) :
    (residueCardAndDegreeObject RD).RecoverableFromAbsGal := by
  sorry

def residueCard_and_degree_recoverable.src : Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.2", sectionId := "prop-1-2" }

/-- 原文の証明文から抽出した、証明が要求するもの(G6)。★下界。 -/
def residueCard_and_degree_recoverable.needs : List ProofObligation :=
  [ .citation "[3] Serre, Local Class Field Theory (Cassels-Frohlich, 1967)"
      "相互律による Γ_K^ab ≅ (K^×)^∧"
      (.inProgress "kbuzzard/ClassFieldTheory"
        "抽象的な reciprocityIso はあるが具体形は未完成。LocalCFT/ は2ファイルのみ") 3,
    .folklore "「it is well-known that」0 → U_K → (K^×)^∧ → Ẑ → 0(典拠なし)" 3,
    .citation "p進対数" "U_K の開部分群を(捩れを除いて)K の開部分群と同一視する"
      (.absent <|
        "mathlib v4.31.0-rc2 全体を再測定(2026-08-14): " ++
        "grep -rliE 'padicLog|Padic.log|padic_log|logOneAdd|log_one_add' → " ++
        "6 件すべて実/複素の解析(Analysis/SpecialFunctions/*・Normed/Module/MultipliableUniformlyOn)。" ++
        "Mathlib/NumberTheory/Padics/ 直下 13 ファイルに log 系の定義は無く、" ++
        "'log' の出現は全て Nat.log / Real.log(付値の対数)。" ++
        "'padicExp|expPadic' → 0 件。" ++
        "公開プロジェクト(ClassFieldTheory・LocalClassFieldTheory)は 2026-08-14 に測定して 0 件" ++
        "(ResearchPaper/lean-ecosystem.json。両リポジトリは本リポジトリに clone されていないので今日は再現していない)。" ++
        "★2026-09-04: 級数の収束のみ `Found/PGC/PadicLog.lean` に自前で構築した(sorry 無し、" ++
        "`padicLog`/`logTerm_summable`)——準同型性(log((1+x)(1+y))=log(1+x)+log(1+y))は" ++
        "まだ無い。★同日、双対の `padicExp`(`Found/PGC/PadicExp.lean`)は加法性" ++
        "`exp(x+y)=exp(x)·exp(y)` まで sorry 無しで到達した(Cauchy 積+二項定理、" ++
        "微分不要)——`exp`/`log` の互逆性が示せれば log の加法性もここから従うが、" ++
        "互逆性自体は未着手。" ++
        "★★★2026-09-04 解消: 在庫確認(`node tools/decl-index.mjs` 再生成後の grep)で、" ++
        "この準同型性が**ℚ_p の場合には既に本リポジトリの別セクションに存在した**と判明" ++
        "(`Found/IUTchIII/PadicLogMul.lean::logOneAdd_mul`——AbsTopIII の log-shell 用に" ++
        "以前構築済み、pGC 側からは未参照だった)。その証明を一般の K へ一般化し、" ++
        "`Found/PGC/PadicLogMul.lean::padicLog_mul` として sorry 無しで完成させた" ++
        "(`K.carrier` の `CharZero`/`IsAddTorsionFree` scoped instance を2つ追加しただけで、" ++
        "証明の筋——形式冪級数 `ABC3.Found.IUTchIII.logOf_mul` の一般形+係数総和の" ++
        "評価橋——はℚ_pの場合と一字一句並行だった)。" ++
        "**この境界外入力(p進対数)は解消——mathlib 不在ではなく、本リポジトリで構築済み。**" ++
        "ただし Proposition 1.2 自体はこれだけでは閉じない: 相互律の同型" ++
        "`Γ_K^ab ≅ (K^×)^∧` そのものと、上の完全系列の典拠は依然として absent") 3 ]

end ABC3.Skeleton.PGC
