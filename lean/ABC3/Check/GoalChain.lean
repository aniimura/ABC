import ABC3

/-!
# ★最終目標からの依存閉包を測る計測器

`ABC3.Skeleton.Goal.abcConjecture_holds` から定数依存グラフを辿り、
`.cache/goal-chain.json` に書き出す。読むのは `tools/goal-chain.mjs`。

★このファイルは `ABC3.lean`(集約)から **import しない**。
全体を import するので、通常のビルドに載せると遅くなるからである。
実行は `node tools/goal-chain.mjs`(内部で `leanfile.mjs` を呼ぶ)。

測るのは 4 つ。

1. **到達** —— 目標から辿れる `ABC3.*` 宣言。
2. **不足** —— そのうち `sorryAx` を直接使うもの(=目標が現在載っている穴)。
3. **過** —— `ABC3.Skeleton.*` のうち目標から**辿れない**もの。
4. **空虚の疑い** —— `ABC3.Interface.*` の構造体のうち、木の中に**住人が 1 つも無い**もの。
   `∀ D : Data, …` 型の主張は、`Data` が空なら自明に真になりうる。
-/

open Lean

namespace ABC3.Check.GoalChain

/-- 定数の本体。★`ConstantInfo.value?` はこの Lean 版では **定理に対して `none`** を返す
(2026-09-09 実測、`lean-idioms.md` #353)。構成子を直接見る。 -/
def valueOf : ConstantInfo -> Option Expr
  | .defnInfo v => some v.value
  | .thmInfo v => some v.value
  | .opaqueInfo v => some v.value
  | _ => none

/-- 定数 `n` が直接使う定数(型と値の両方から集める)。 -/
def directDeps (env : Environment) (n : Name) : Array Name :=
  match env.find? n with
  | none => #[]
  | some ci =>
    match valueOf ci with
    | some v => ci.type.getUsedConstants ++ v.getUsedConstants
    | none => ci.type.getUsedConstants

/-- 推移閉包。 -/
partial def closureAux (env : Environment) : List Name → NameSet → NameSet
  | [], acc => acc
  | n :: rest, acc =>
    if acc.contains n then closureAux env rest acc
    else closureAux env ((directDeps env n).toList ++ rest) (acc.insert n)

/-- 前置きが `p` で始まるか。 -/
def hasPrefix (n : Name) (p : Name) : Bool := p.isPrefixOf n

/-- `∀`-束縛を剥がした結論の先頭定数。 -/
partial def conclHead : Expr → Option Name
  | .forallE _ _ b _ => conclHead b
  | e => match e.getAppFn with
         | .const c _ => some c
         | _ => none

/-- 主張として数えない接尾辞(出典メタデータ・検査用の欄)。 -/
def isMeta (s : String) : Bool :=
  s.endsWith ".src" || s.endsWith ".needs" || s.endsWith ".inventory"
    || s.endsWith ".nonvacuous" || s.endsWith ".negControl"
    || s.endsWith ".loadBearing" || s.endsWith ".waiting"
    || s.endsWith ".deviation" || s.endsWith ".note"

/-- ★Lean が自動生成した宣言か。手で書いた主張ではないので孤児の数から除く。

実測(2026-09-09): 孤児 156 件のうち **33 件**がこれ——`.mk.injEq` / `.recOn` /
`.ctorIdx` / `.sizeOf_spec` / `.eq_1`、および `structure` の射影
(`FilteredGroup.antitone` など)。構造体を消さない限り消せないので、
「消費者がいない」の対象ではない。 -/
def isGenerated (env : Environment) (n : Name) : Bool :=
  let s := n.toString
  let gen : List String :=
    [".injEq", ".sizeOf_spec", ".ctorIdx", ".recOn", ".rec", ".casesOn",
     ".noConfusion", ".noConfusionType", ".below", ".brecOn", ".ndrec",
     ".mk", ".eq_def", ".sizeOf_inst", ".ofNat", ".toCtorIdx"]
  if gen.any (fun t => s.endsWith t) then true
  else
    -- `foo.eq_1`, `foo.eq_2`, … (等式補題)
    let last := n.getString!
    if last.startsWith "eq_" && (last.drop 3).all Char.isDigit && last.length > 3 then true
    else
      -- `structure` の射影
      match n with
      | .str par f => isStructure env par && (getStructureFields env par).contains (Name.mkSimple f)
      | _ => false

/-- 型に `∀`/`→` が無い(＝閉じた住人)か。 -/
def isClosedType : Expr -> Bool
  | .forallE .. => false
  | _ => true

/-- 計測本体。JSON 文字列を返す。 -/
def report (env : Environment) : String := Id.run do
  let goal : Name := `ABC3.Skeleton.Goal.abcConjecture_holds
  let cl := closureAux env [goal] {}
  -- ① 到達した ABC3 宣言 / ② そのうち sorry を直接使うもの
  let mut reached : Array Name := #[]
  let mut sorries : Array Name := #[]
  for n in cl.toList do
    if hasPrefix n `ABC3 && !n.isInternal then
      reached := reached.push n
      if (directDeps env n).contains ``sorryAx then
        sorries := sorries.push n
  -- 全 ABC3 定数を 1 度だけ走査する
  let mut skelAll : Array Name := #[]
  let mut used : NameSet := {}
  let mut inhabAny : Std.HashMap Name Nat := {}
  let mut inhabClosed : Std.HashMap Name (Array Name) := {}
  let mut abc3Count : Nat := 0
  for (n, ci) in env.constants.toList do
    if n.isInternal then continue
    if !(hasPrefix n `ABC3) then continue
    abc3Count := abc3Count + 1
    for d in directDeps env n do
      if d != n then used := used.insert d
    if hasPrefix n `ABC3.Skeleton && !(isMeta n.toString) && !(isGenerated env n) then
      skelAll := skelAll.push n
    match conclHead ci.type with
    | some h =>
      if hasPrefix h `ABC3.Interface && isStructure env h && !(h.isPrefixOf n) then
        inhabAny := inhabAny.insert h ((inhabAny.getD h 0) + 1)
        if isClosedType ci.type then
          inhabClosed := inhabClosed.insert h ((inhabClosed.getD h #[]).push n)
    | none => pure ()
  -- ③ Interface の構造体一覧
  let mut ifaces : Array (Name × Nat × Array Name) := #[]
  for (n, _) in env.constants.toList do
    if n.isInternal then continue
    if hasPrefix n `ABC3.Interface && isStructure env n then
      ifaces := ifaces.push (n, inhabAny.getD n 0, inhabClosed.getD n #[])
  let unreached := skelAll.filter (fun n => !cl.contains n)
  -- ④ 孤児 —— ABC3 のどの宣言からも参照されていない Skeleton の主張
  --    (★目標の sorry がどこで切れているかに依存しない指標)
  let orphans := skelAll.filter (fun n => !used.contains n)
  let q := fun (n : Name) => "\"" ++ n.toString ++ "\""
  let arr := fun (a : Array Name) =>
    "[" ++ String.intercalate "," (a.qsort (·.toString < ·.toString) |>.toList.map q) ++ "]"
  let tri := fun (x : Name × Nat × Array Name) =>
    "{\"name\":" ++ q x.1 ++ ",\"inhabitants\":" ++ toString x.2.1
      ++ ",\"closed\":" ++ toString x.2.2.size
      ++ ",\"closedNames\":" ++ arr x.2.2 ++ "}"
  return "{
  \"goal\": " ++ q goal
    ++ ",
  \"abc3Constants\": " ++ toString abc3Count
    ++ ",
  \"reached\": " ++ arr reached
    ++ ",
  \"sorryLeaves\": " ++ arr sorries
    ++ ",
  \"skeletonAll\": " ++ arr skelAll
    ++ ",
  \"skeletonUnreached\": " ++ arr unreached
    ++ ",
  \"skeletonOrphans\": " ++ arr orphans
    ++ ",
  \"interfaces\": [" ++ String.intercalate ","
         (ifaces.qsort (fun a b => a.1.toString < b.1.toString) |>.toList.map tri) ++ "]
}
"

end ABC3.Check.GoalChain

open Lean in
#eval show CoreM Unit from do
  let env ← getEnv
  let s := ABC3.Check.GoalChain.report env
  IO.FS.createDirAll ".cache"
  IO.FS.writeFile ".cache/goal-chain.json" s
  IO.println s!"goal-chain: wrote .cache/goal-chain.json ({s.length} bytes)"
